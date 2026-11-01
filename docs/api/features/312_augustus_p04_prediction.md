# AUGUSTUS：基因预测、训练与流水线

`dnakit.annotation` 中的这组函数是 AUGUSTUS 3.5.0 程序与辅助脚本（[AUGUSTUS](https://github.com/Gaius-Augustus/Augustus)）的类型化封装。每个函数用显式参数列表运行一个固定程序，不经过 shell。请先安装 AUGUSTUS（`conda install -c bioconda augustus`）并把其 `bin` 目录放入 PATH，或传入 `augustus_bin_dir`；缺少程序时抛出带安装提示的 `BackendUnavailableError`。物种参数集从 `AUGUSTUS_CONFIG_PATH`（或 `config_path`）读取，训练类函数会写入该目录，因此需可写。输出路径必须是新路径，不会覆盖已有文件。每个结果是包含 `status`、`backend`、`version`、`command`、`return_code`、`outputs`、`log_tail`、`provenance` 以及少量解析计数的字典。

## 1) `dnakit.annotation.augustus_auto_predict()`

- **作用：** 分块的全基因组基因预测并合并结果（autoAugPred.pl --singleCPU）。
- **处理方法：** 通过 DNAKit 受限执行器运行`augustus_options`、`cpus`、`extrinsic_config_path`、`hints_path`，使用固定参数白名单。
- **API：** `dnakit.annotation.augustus_auto_predict(genome_path, species, work_dir, *, hints_path=None, extrinsic_config_path=None, utr=False, augustus_options=(), use_existing=False, cpus=1, config_path=None, augustus_bin_dir=None, verbose=True)`。
- **输入：** 必需参数：`genome_path`、`species`、`work_dir`；可选关键字参数：`hints_path`、`extrinsic_config_path`、`utr`、`augustus_options`、`use_existing`、`cpus`、`config_path`、`augustus_bin_dir`。输入文件必须存在，输出路径必须是新路径。
- **输出：** 含 `status`、`backend`、`version`、`command`、`return_code`、`outputs`、`log_tail`、`provenance` 以及 `lines`、`records`、`genes`、`transcripts` 等解析计数的字典。
- **示例代码：** 此处从略；示例需要体积过大的数据，或 AUGUSTUS 软件包之外的外部工具（见测试与 manifest 说明）。

## 2) `dnakit.annotation.augustus_auto_train_and_predict()`

- **作用：** 自动训练并注释的流水线（autoAug.pl --singleCPU）。
- **处理方法：** 通过 DNAKit 受限执行器运行`cdna_path`、`cpus`、`est_alignment_path`、`genome_path`，使用固定参数白名单。
- **API：** `dnakit.annotation.augustus_auto_train_and_predict(genome_path, species, work_dir, *, training_path=None, hints_path=None, est_alignment_path=None, cdna_path=None, no_utr=False, optimization_rounds=None, max_intron_length=None, use_existing=False, cpus=1, config_path=None, augustus_bin_dir=None, verbose=True)`。
- **输入：** 必需参数：`genome_path`、`species`、`work_dir`；可选关键字参数：`training_path`、`hints_path`、`est_alignment_path`、`cdna_path`、`no_utr`、`optimization_rounds`、`max_intron_length`、`use_existing`、`cpus`、`config_path`、`augustus_bin_dir`。输入文件必须存在，输出路径必须是新路径。
- **输出：** 含 `status`、`backend`、`version`、`command`、`return_code`、`outputs`、`log_tail`、`provenance` 以及 `lines`、`records`、`genes`、`transcripts` 等解析计数的字典。
- **示例代码：** 此处从略；示例需要体积过大的数据，或 AUGUSTUS 软件包之外的外部工具（见测试与 manifest 说明）。

## 3) `dnakit.annotation.augustus_choose_parameter_set()`

- **作用：** 为每条序列选择 AUGUSTUS 参数集并预测基因（augustify.py）。
- **处理方法：** 通过 DNAKit 受限执行器运行`augustify.py`、`classification_path`、`prediction_path`、`species`，使用固定参数白名单。
- **API：** `dnakit.annotation.augustus_choose_parameter_set(genome_path, species_list, *, mode='dominant_species', classification_path=None, prediction_path=None, threads=1, config_path=None, augustus_bin_dir=None, verbose=True)`。
- **输入：** 必需参数：`genome_path`、`species_list`；可选关键字参数：`mode`、`classification_path`、`prediction_path`、`threads`、`config_path`、`augustus_bin_dir`。输入文件必须存在，输出路径必须是新路径。
- **输出：** 含 `status`、`backend`、`version`、`command`、`return_code`、`outputs`、`log_tail`、`provenance` 以及 `lines`、`records`、`genes`、`transcripts` 等解析计数的字典。
- **示例代码：** 此处从略；示例需要体积过大的数据，或 AUGUSTUS 软件包之外的外部工具（见测试与 manifest 说明）。

## 4) `dnakit.annotation.augustus_fix_in_frame_stop_genes()`

- **作用：** 用 MEA 重新预测含框内终止密码子的基因并替换（fix_in_frame_stop_codon_genes.py）。
- **处理方法：** 通过 DNAKit 受限执行器运行`.gff3`、`additional_augustus_arguments`、`annotation_path`、`bad_genes_path`，使用固定参数白名单。
- **API：** `dnakit.annotation.augustus_fix_in_frame_stop_genes(genome_path, annotation_path, bad_genes_path, output_stem, species, *, annotation_format='gtf', hints_path=None, extrinsic_config_path=None, softmasking=False, utr=False, print_utr=False, additional_augustus_arguments=None, cores=1, config_path=None, augustus_bin_dir=None, verbose=True)`。
- **输入：** 必需参数：`genome_path`、`annotation_path`、`bad_genes_path`、`output_stem`、`species`；可选关键字参数：`annotation_format`、`hints_path`、`extrinsic_config_path`、`softmasking`、`utr`、`print_utr`、`additional_augustus_arguments`、`cores`、`config_path`、`augustus_bin_dir`。输入文件必须存在，输出路径必须是新路径。
- **输出：** 含 `status`、`backend`、`version`、`command`、`return_code`、`outputs`、`log_tail`、`provenance` 以及 `lines`、`records`、`genes`、`transcripts` 等解析计数的字典。
- **示例代码：** 此处从略；示例需要体积过大的数据，或 AUGUSTUS 软件包之外的外部工具（见测试与 manifest 说明）。

## 5) `dnakit.annotation.augustus_homologous_gene_mapping()`

- **作用：** 用 HAL 比对在不同基因组间映射基因模型（homGeneMapping）。
- **处理方法：** 通过 DNAKit 受限执行器运行`gene_files_table_path`、`gtffilenames.tbl`、`halLiftover`、`hal_bin_dir`，使用固定参数白名单。
- **API：** `dnakit.annotation.augustus_homologous_gene_mapping(hal_path, gene_files_table_path, output_dir, *, no_dupes=False, details=False, unmapped=False, homologs_path=None, hal_bin_dir=None, cpus=1, augustus_bin_dir=None, verbose=True)`。
- **输入：** 必需参数：`hal_path`、`gene_files_table_path`、`output_dir`；可选关键字参数：`no_dupes`、`details`、`unmapped`、`homologs_path`、`hal_bin_dir`、`cpus`、`augustus_bin_dir`。输入文件必须存在，输出路径必须是新路径。
- **输出：** 含 `status`、`backend`、`version`、`command`、`return_code`、`outputs`、`log_tail`、`provenance` 以及 `lines`、`records`、`genes`、`transcripts` 等解析计数的字典。
- **示例代码：** 此处从略；示例需要体积过大的数据，或 AUGUSTUS 软件包之外的外部工具（见测试与 manifest 说明）。

## 6) `dnakit.annotation.augustus_join_gene_sets()`

- **作用：** 合并多个 GTF 基因集并选择最佳转录本（joingenes）。
- **处理方法：** 通过 DNAKit 受限执行器运行`alternatives`、`bacterium`、`error_distance`、`eukaryote`，使用固定参数白名单。
- **API：** `dnakit.annotation.augustus_join_gene_sets(gene_set_paths, output_path, *, priorities=None, error_distance=None, gene_model=None, alternatives=False, suppress_priorities=None, stop_in_coding=False, no_join=False, no_selection=False, augustus_bin_dir=None, verbose=True)`。
- **输入：** 必需参数：`gene_set_paths`、`output_path`；可选关键字参数：`priorities`、`error_distance`、`gene_model`、`alternatives`、`suppress_priorities`、`stop_in_coding`、`no_join`、`no_selection`、`augustus_bin_dir`。输入文件必须存在，输出路径必须是新路径。
- **输出：** 含 `status`、`backend`、`version`、`command`、`return_code`、`outputs`、`log_tail`、`provenance` 以及 `lines`、`records`、`genes`、`transcripts` 等解析计数的字典。
- **示例代码：**

```python
import random, tempfile
from pathlib import Path
from dnakit.annotation import augustus_p04_prediction_functions as m
w = Path(tempfile.mkdtemp())
(w / "genes.gtf").write_text(
    'chr1\tAUGUSTUS\tgene\t101\t400\t.\t+\t.\tg1\n'
    'chr1\tAUGUSTUS\ttranscript\t101\t400\t.\t+\t.\tg1.t1\n'
    'chr1\tAUGUSTUS\tstart_codon\t101\t103\t.\t+\t0\ttranscript_id "g1.t1"; gene_id "g1";\n'
    'chr1\tAUGUSTUS\tCDS\t101\t200\t.\t+\t0\ttranscript_id "g1.t1"; gene_id "g1";\n'
    'chr1\tAUGUSTUS\tCDS\t301\t400\t.\t+\t1\ttranscript_id "g1.t1"; gene_id "g1";\n'
)
r = m.augustus_join_gene_sets([str(w / "genes.gtf")], str(w / "merged.gtf"), verbose=False)
print(r["genes"], r["transcripts"])
```

- **示例输出：**

```text
1 1
```

## 7) `dnakit.annotation.augustus_join_segment_predictions()`

- **作用：** 合并重叠序列片段的 AUGUSTUS 预测结果（join_aug_pred.pl）。
- **处理方法：** 通过 DNAKit 受限执行器运行`drop_list_path`、`join_aug_pred.pl`、`prediction_paths`，使用固定参数白名单。
- **API：** `dnakit.annotation.augustus_join_segment_predictions(prediction_paths, output_path, *, drop_list_path=None, augustus_bin_dir=None, verbose=True)`。
- **输入：** 必需参数：`prediction_paths`、`output_path`；可选关键字参数：`drop_list_path`、`augustus_bin_dir`。输入文件必须存在，输出路径必须是新路径。
- **输出：** 含 `status`、`backend`、`version`、`command`、`return_code`、`outputs`、`log_tail`、`provenance` 以及 `lines`、`records`、`genes`、`transcripts` 等解析计数的字典。
- **示例代码：**

```python
import random, tempfile
from pathlib import Path
from dnakit.annotation import augustus_p04_prediction_functions as m
w = Path(tempfile.mkdtemp())
rng = random.Random(7)
(w / "genome.fa").write_text(">chr1\n" + "".join(rng.choices("ACGT", k=4000)) + "\n")
a = m.augustus_predict_genes(str(w / "genome.fa"), "human", str(w / "a.gtf"), verbose=False)
j = m.augustus_join_segment_predictions([a["outputs"]["predictions"]] * 2, str(w / "joined.gtf"), verbose=False)
print(j["input_files"], j["genes"] == a["genes"])
```

- **示例输出：**

```text
2 True
```

## 8) `dnakit.annotation.augustus_optimize_parameters()`

- **作用：** 用 k 折交叉验证自动优化 AUGUSTUS 元参数（optimize_augustus.pl）。
- **处理方法：** 通过 DNAKit 受限执行器运行`augustus_train_species`、`only_train_path`、`optimize_augustus.pl`、`species`，使用固定参数白名单。
- **API：** `dnakit.annotation.augustus_optimize_parameters(species, training_path, work_dir, *, rounds=5, cpus=1, kfold=8, pstep=5, metapars_path=None, only_train_path=None, utr=False, sens_spec_bias=None, train_only_utr=False, no_train_pars=False, cleanup=False, config_path=None, augustus_bin_dir=None, verbose=True)`。
- **输入：** 必需参数：`species`、`training_path`、`work_dir`；可选关键字参数：`rounds`、`cpus`、`kfold`、`pstep`、`metapars_path`、`only_train_path`、`utr`、`sens_spec_bias`、`train_only_utr`、`no_train_pars`、`cleanup`、`config_path`、`augustus_bin_dir`。输入文件必须存在，输出路径必须是新路径。
- **输出：** 含 `status`、`backend`、`version`、`command`、`return_code`、`outputs`、`log_tail`、`provenance` 以及 `lines`、`records`、`genes`、`transcripts` 等解析计数的字典。
- **示例代码：** 此处从略；示例需要体积过大的数据，或 AUGUSTUS 软件包之外的外部工具（见测试与 manifest 说明）。

## 9) `dnakit.annotation.augustus_predict_comparative()`

- **作用：** 在多物种基因组比对（MAF + 树）上做比较基因预测（CGP）。
- **处理方法：** 通过 DNAKit 受限执行器运行`alignment_path`、`dd_iterations`、`dual_decomposition`、`output_path`，使用固定参数白名单。
- **API：** `dnakit.annotation.augustus_predict_comparative(species, tree_path, alignment_path, species_filenames_path, output_path, *, dual_decomposition=False, dd_iterations=None, only_species=None, softmasking=None, config_path=None, augustus_bin_dir=None, verbose=True)`。
- **输入：** 必需参数：`species`、`tree_path`、`alignment_path`、`species_filenames_path`、`output_path`；可选关键字参数：`dual_decomposition`、`dd_iterations`、`only_species`、`softmasking`、`config_path`、`augustus_bin_dir`。输入文件必须存在，输出路径必须是新路径。
- **输出：** 含 `status`、`backend`、`version`、`command`、`return_code`、`outputs`、`log_tail`、`provenance` 以及 `lines`、`records`、`genes`、`transcripts` 等解析计数的字典。
- **示例代码：** 此处从略；示例需要体积过大的数据，或 AUGUSTUS 软件包之外的外部工具（见测试与 manifest 说明）。

## 10) `dnakit.annotation.augustus_predict_genes()`

- **作用：** 用 AUGUSTUS 在基因组 DNA 上预测基因结构（HMM，可选 hints；含抽样候选转录本、MEA 解码、蛋白家族 profile 引导三种模式）。
- **处理方法：** 通过 DNAKit 受限执行器运行`augustus_extract_annotation_fasta`、`coding_sequences`、`extrinsic_config_path`、`hints_path`，使用固定参数白名单。
- **API：** `dnakit.annotation.augustus_predict_genes(sequence_path, species, output_path, *, strand='both', gene_model='partial', single_strand=False, hints_path=None, extrinsic_config_path=None, utr=False, softmasking=None, alternatives_from_evidence=None, alternatives_from_sampling=False, sample=None, min_exon_intron_prob=None, min_mean_exon_intron_prob=None, max_tracks=None, mea=False, protein_profile_path=None, gff3=False, protein_sequences=False, coding_sequences=False, prediction_start=None, prediction_end=None, no_in_frame_stop=False, unique_gene_id=False, config_path=None, augustus_bin_dir=None, verbose=True)`。
- **输入：** 必需参数：`sequence_path`、`species`、`output_path`；可选关键字参数：`strand`、`gene_model`、`single_strand`、`hints_path`、`extrinsic_config_path`、`utr`、`softmasking`、`alternatives_from_evidence`、`alternatives_from_sampling`、`sample`、`min_exon_intron_prob`、`min_mean_exon_intron_prob`、`max_tracks`、`mea`、`protein_profile_path`、`gff3`、`protein_sequences`、`coding_sequences`、`prediction_start`、`prediction_end`、`no_in_frame_stop`、`unique_gene_id`、`config_path`、`augustus_bin_dir`。输入文件必须存在，输出路径必须是新路径。
- **输出：** 含 `status`、`backend`、`version`、`command`、`return_code`、`outputs`、`log_tail`、`provenance` 以及 `lines`、`records`、`genes`、`transcripts` 等解析计数的字典。
- **示例代码：**

```python
import random, tempfile
from pathlib import Path
from dnakit.annotation import augustus_p04_prediction_functions as m
w = Path(tempfile.mkdtemp())
rng = random.Random(7)
(w / "genome.fa").write_text(">chr1\n" + "".join(rng.choices("ACGT", k=4000)) + "\n")
r = m.augustus_predict_genes(str(w / "genome.fa"), "human", str(w / "pred.gtf"), mea=True, verbose=False)
print(r["status"], r["version"], r["genes"] >= 0)
print(r["command"][1:5])
```

- **示例输出：**

```text
completed 3.5.0 True
['--species=human', '--strand=both', '--genemodel=partial', '--UTR=off']
```

## 11) `dnakit.annotation.augustus_predict_utrs_from_rnaseq()`

- **作用：** 用 RNA-Seq 覆盖度预测已注释基因的 UTR（utrrnaseq）。
- **处理方法：** 通过 DNAKit 受限执行器运行`coding_regions_path`、`introns_path`、`output_path`、`repeats_path`，使用固定参数白名单。
- **API：** `dnakit.annotation.augustus_predict_utrs_from_rnaseq(genome_path, coding_regions_path, introns_path, wiggle_path, output_path, *, repeats_path=None, smoothing_window_size=150, read_length=150, limit=5000, drop_window_size=50, minimum_length=2, minimum_average_coverage=10, percent_window=0.6, percent_intron=0.5, percent_multiplicity=0.1, splice_sites='GT_AG', zero_coverage=False, augustus_bin_dir=None, verbose=True)`。
- **输入：** 必需参数：`genome_path`、`coding_regions_path`、`introns_path`、`wiggle_path`、`output_path`；可选关键字参数：`repeats_path`、`smoothing_window_size`、`read_length`、`limit`、`drop_window_size`、`minimum_length`、`minimum_average_coverage`、`percent_window`、`percent_intron`、`percent_multiplicity`、`splice_sites`、`zero_coverage`、`augustus_bin_dir`。输入文件必须存在，输出路径必须是新路径。
- **输出：** 含 `status`、`backend`、`version`、`command`、`return_code`、`outputs`、`log_tail`、`provenance` 以及 `lines`、`records`、`genes`、`transcripts` 等解析计数的字典。
- **示例代码：** 此处从略；示例需要体积过大的数据，或 AUGUSTUS 软件包之外的外部工具（见测试与 manifest 说明）。

## 12) `dnakit.annotation.augustus_protein_profile_block_search()`

- **作用：** 用 fastBlockSearch 在基因组 DNA 中定位蛋白家族 profile 的 block 命中。
- **处理方法：** 通过 DNAKit 受限执行器运行`.prfl`、`fastBlockSearch`、`msa2prfl.pl`、`output_path`，使用固定参数白名单。
- **API：** `dnakit.annotation.augustus_protein_profile_block_search(sequence_path, profile_path, species, output_path, *, config_path=None, augustus_bin_dir=None, verbose=True)`。
- **输入：** 必需参数：`sequence_path`、`profile_path`、`species`、`output_path`；可选关键字参数：`config_path`、`augustus_bin_dir`。输入文件必须存在，输出路径必须是新路径。
- **输出：** 含 `status`、`backend`、`version`、`command`、`return_code`、`outputs`、`log_tail`、`provenance` 以及 `lines`、`records`、`genes`、`transcripts` 等解析计数的字典。
- **示例代码：** 此处从略；示例需要体积过大的数据，或 AUGUSTUS 软件包之外的外部工具（见测试与 manifest 说明）。

## 13) `dnakit.annotation.augustus_train_species()`

- **作用：** 由训练基因训练新物种的基因模型参数（etraining 或 autoAugTrain.pl）。
- **处理方法：** 通过 DNAKit 受限执行器运行`autoAugTrain.pl`、`cpus`、`crf`、`est_alignment_path`，使用固定参数白名单。
- **API：** `dnakit.annotation.augustus_train_species(species, training_path, *, method='etraining', prokaryotic=False, set_stop_codon_frequencies=False, genome_path=None, utr=False, est_alignment_path=None, work_dir=None, optimization_rounds=None, crf=False, flanking_dna=None, use_existing=False, cpus=None, config_path=None, augustus_bin_dir=None, verbose=True)`。
- **输入：** 必需参数：`species`、`training_path`；可选关键字参数：`method`、`prokaryotic`、`set_stop_codon_frequencies`、`genome_path`、`utr`、`est_alignment_path`、`work_dir`、`optimization_rounds`、`crf`、`flanking_dna`、`use_existing`、`cpus`、`config_path`、`augustus_bin_dir`。输入文件必须存在，输出路径必须是新路径。
- **输出：** 含 `status`、`backend`、`version`、`command`、`return_code`、`outputs`、`log_tail`、`provenance` 以及 `lines`、`records`、`genes`、`transcripts` 等解析计数的字典。
- **示例代码：** 此处从略；示例需要体积过大的数据，或 AUGUSTUS 软件包之外的外部工具（见测试与 manifest 说明）。

