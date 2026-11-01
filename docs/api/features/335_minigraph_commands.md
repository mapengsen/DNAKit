# minigraph：序列到图的映射、图构建、气泡调用与结构变异脚本

`dnakit.alignment` 中的这一组函数是对 [minigraph](https://github.com/lh3/minigraph) 及其配套脚本的命令行封装。`map_sequences_to_graph`、`graph_segment_coverage`、`build_sequence_graph` 和 `call_graph_bubbles` 运行 `minigraph` 程序（映射、`--cov`、`-cxggs`、`--call`）。其余函数通过 `k8` JavaScript shell 运行 `mgutils.js`（重命名 FASTA、paf2bl、anno、anno2tbl、stableGaf、extractseg、merge、segfreq、genecopy）、`mgutils-es6.js`（merge2vcf、addsample）或 `gafcall.js`（extract、merge、eval、view、join）的一个子命令。未提供 `mgutils.js joinfa`：它在所有已发布的 k8 上都会失败（脚本用 `Bytes` 参数调用 `Bytes.set`，k8 1.2 不接受）。

测试版本为 minigraph 0.21 和 k8 1.2（MIT 许可，作为独立程序启动）。用 `conda install -c bioconda minigraph k8` 安装或按项目页面编译；把 `minigraph`、`k8` 和脚本（源码 `misc/` 目录中的 `mgutils.js`、`mgutils-es6.js`、`gafcall.js`）放到 PATH 上，或传入 `minigraph_path`、`k8_path`、`mgutils_path`、`mgutils_es6_path`、`gafcall_path`。每个函数都写入新的 `output_directory`；输出到标准输出的程序会被写入其中的文件，标准错误写入 `tool.log`。值为 `None` 的选项不会传出，因此使用工具默认值。示例使用 `verbose=False`。

## 示例数据与示例

**示例代码：**

```python
from pathlib import Path
import random, tempfile
from dnakit.alignment import minigraph_p11_cli_functions as mg

rng = random.Random(3)
work = Path(tempfile.mkdtemp())
ref = "".join(rng.choice("ACGT") for _ in range(8000))
alt = ref[:3000] + "".join(rng.choice("ACGT") for _ in range(400)) + ref[3000:]
(work / "ref.fa").write_text(f">ref\n{ref}\n")
(work / "alt.fa").write_text(f">alt\n{alt}\n")
graph = mg.build_sequence_graph(work / "ref.fa", [work / "alt.fa"], work / "graph", min_aln_length=1000, min_depth_aln_length=1000, threads=1, verbose=False)
print(graph["version"], graph["segments"], graph["links"])
calls = mg.call_graph_bubbles(graph["outputs"]["output"], work / "alt.fa", work / "calls", min_aln_length=1000, min_depth_aln_length=1000, threads=1, verbose=False)
print(calls["bubbles"], calls["rows"][0][:3])
mapped = mg.map_sequences_to_graph(graph["outputs"]["output"], [work / "alt.fa"], work / "map", preset="asm", base_alignment=True, threads=1, verbose=False)
print(mapped["records"], mapped["mapped_queries"])
```

**示例输出：**

```text
0.21 3 3
1 ['ref', '2999', '2999']
1 1
```

**示例代码：**

```python
from pathlib import Path
import tempfile
from dnakit.alignment import minigraph_p11_cli_functions as mg

work = Path(tempfile.mkdtemp())
(work / "ins.bed").write_text("chr1\t100\t300\t200\nchr1\t500\t900\t400\n")
(work / "ins.paf").write_text("q1\t800\t0\t800\t+\tchr1\t10000\t100\t900\t780\t800\t60\ttp:A:P\tde:f:0.05\n")
print(mg.blacklist_insert_alignments(work / "ins.paf", work / "bl", verbose=False)["regions"])
anno = mg.annotate_inserted_sequences(work / "ins.bed", work / "anno", verbose=False)
print(anno["annotated"], mg.summarize_insert_annotation(anno["outputs"]["output"], work / "tbl", verbose=False)["classes"][0]["label"])
```

**示例输出：**

```text
1
2 Other-repeat
```


## 1) `dnakit.alignment.minigraph_p11_cli_functions.add_sample_names_to_calls()`

- **作用：** 给合并的气泡结果表头加入样本名（mgutils-es6.js addsample）。
- **处理方法：** 运行固定的参数列表（不经过 shell）；参数含义、取值范围和工具默认值见英文版同名小节的 Method 说明。
- **API：** `dnakit.alignment.minigraph_p11_cli_functions.add_sample_names_to_calls(merged_calls: 'str | Path', samples: 'Sequence[str]', output_directory: 'str | Path', *, output_name: 'str' = 'merged_named.txt', mgutils_es6_path: 'str | None' = None, k8_path: 'str | None' = None, timeout_seconds: 'float' = 3600.0, verbose: 'bool' = True) -> 'dict[str, Any]'`。
- **输入：** 必填：`merged_calls`, `samples`, `output_directory`；其余参数均为仅关键字参数且可选（`None` 或 `False` 表示使用工具自己的默认值）。
- **输出：** 含 `status, backend, operation, version, command, output_directory, outputs, log_tail, provenance` 的字典，另有该函数特有的键（见英文版 Output）。

## 2) `dnakit.alignment.minigraph_p11_cli_functions.annotate_inserted_sequences()`

- **作用：** 注释较短的插入序列（重复、低复杂度、缺口、片段重复；mgutils.js anno）。
- **处理方法：** 运行固定的参数列表（不经过 shell）；参数含义、取值范围和工具默认值见英文版同名小节的 Method 说明。
- **API：** `dnakit.alignment.minigraph_p11_cli_functions.annotate_inserted_sequences(bed: 'str | Path', output_directory: 'str | Path', *, output_name: 'str' = 'anno.tsv', repeatmasker_out: 'str | Path | None' = None, gaps: 'str | Path | None' = None, low_complexity: 'str | Path | None' = None, tandem_repeats: 'str | Path | None' = None, paf: 'str | Path | None' = None, centromere: 'str | Path | None' = None, bubbles: 'str | Path | None' = None, segdup: 'str | Path | None' = None, min_feature_length: 'int | None' = None, min_segdup_fraction: 'float | None' = None, annotate_mobile: 'bool' = False, mgutils_path: 'str | None' = None, k8_path: 'str | None' = None, timeout_seconds: 'float' = 3600.0, verbose: 'bool' = True) -> 'dict[str, Any]'`。
- **输入：** 必填：`bed`, `output_directory`；其余参数均为仅关键字参数且可选（`None` 或 `False` 表示使用工具自己的默认值）。
- **输出：** 含 `status, backend, operation, version, command, output_directory, outputs, log_tail, provenance` 的字典，另有该函数特有的键（见英文版 Output）。

## 3) `dnakit.alignment.minigraph_p11_cli_functions.blacklist_insert_alignments()`

- **作用：** 由插入序列到参考的比对得到要列入黑名单的参考区域（mgutils.js paf2bl）。
- **处理方法：** 运行固定的参数列表（不经过 shell）；参数含义、取值范围和工具默认值见英文版同名小节的 Method 说明。
- **API：** `dnakit.alignment.minigraph_p11_cli_functions.blacklist_insert_alignments(paf: 'str | Path', output_directory: 'str | Path', *, output_name: 'str' = 'blacklist.bed', min_divergence: 'float' = 0.01, identical_only: 'bool' = False, mgutils_path: 'str | None' = None, k8_path: 'str | None' = None, timeout_seconds: 'float' = 3600.0, verbose: 'bool' = True) -> 'dict[str, Any]'`。
- **输入：** 必填：`paf`, `output_directory`；其余参数均为仅关键字参数且可选（`None` 或 `False` 表示使用工具自己的默认值）。
- **输出：** 含 `status, backend, operation, version, command, output_directory, outputs, log_tail, provenance` 的字典，另有该函数特有的键（见英文版 Output）。

## 4) `dnakit.alignment.minigraph_p11_cli_functions.build_sequence_graph()`

- **作用：** 把组装逐个加入参考，增量构建序列图（minigraph -cxggs）。
- **处理方法：** 运行固定的参数列表（不经过 shell）；参数含义、取值范围和工具默认值见英文版同名小节的 Method 说明。
- **API：** `dnakit.alignment.minigraph_p11_cli_functions.build_sequence_graph(reference: 'str | Path', assemblies: 'Sequence[str | Path]', output_directory: 'str | Path', *, output_name: 'str' = 'graph.gfa', preset: 'str | None' = 'ggs', base_alignment: 'bool' = True, min_mapq: 'int | None' = None, min_aln_length: 'int | None' = None, min_depth_aln_length: 'int | None' = None, min_variant_length: 'int | None' = None, kmer_size: 'int | None' = None, window_size: 'int | None' = None, threads: 'int' = 4, minigraph_path: 'str | None' = None, timeout_seconds: 'float' = 3600.0, verbose: 'bool' = True) -> 'dict[str, Any]'`。
- **输入：** 必填：`reference`, `assemblies`, `output_directory`；其余参数均为仅关键字参数且可选（`None` 或 `False` 表示使用工具自己的默认值）。
- **输出：** 含 `status, backend, operation, version, command, output_directory, outputs, log_tail, provenance` 的字典，另有该函数特有的键（见英文版 Output）。

## 5) `dnakit.alignment.minigraph_p11_cli_functions.call_graph_bubbles()`

- **作用：** 在图的每个气泡中给出某组装走过的路径（minigraph --call）。
- **处理方法：** 运行固定的参数列表（不经过 shell）；参数含义、取值范围和工具默认值见英文版同名小节的 Method 说明。
- **API：** `dnakit.alignment.minigraph_p11_cli_functions.call_graph_bubbles(graph: 'str | Path', assembly: 'str | Path', output_directory: 'str | Path', *, output_name: 'str' = 'calls.bed', preset: 'str | None' = 'asm', base_alignment: 'bool' = True, min_mapq: 'int | None' = None, min_aln_length: 'int | None' = None, min_depth_aln_length: 'int | None' = None, min_variant_length: 'int | None' = None, kmer_size: 'int | None' = None, window_size: 'int | None' = None, threads: 'int' = 4, minigraph_path: 'str | None' = None, timeout_seconds: 'float' = 3600.0, verbose: 'bool' = True) -> 'dict[str, Any]'`。
- **输入：** 必填：`graph`, `assembly`, `output_directory`；其余参数均为仅关键字参数且可选（`None` 或 `False` 表示使用工具自己的默认值）。
- **输出：** 含 `status, backend, operation, version, command, output_directory, outputs, log_tail, provenance` 的字典，另有该函数特有的键（见英文版 Output）。

## 6) `dnakit.alignment.minigraph_p11_cli_functions.convert_calls_to_vcf()`

- **作用：** 把合并的气泡结果转换为 VCF（mgutils-es6.js merge2vcf）。
- **处理方法：** 运行固定的参数列表（不经过 shell）；参数含义、取值范围和工具默认值见英文版同名小节的 Method 说明。
- **API：** `dnakit.alignment.minigraph_p11_cli_functions.convert_calls_to_vcf(merged_calls: 'str | Path', output_directory: 'str | Path', *, output_name: 'str' = 'calls.vcf', reference_index: 'int | None' = None, max_alleles: 'int | None' = None, samples: 'Sequence[str] | None' = None, mgutils_es6_path: 'str | None' = None, k8_path: 'str | None' = None, timeout_seconds: 'float' = 3600.0, verbose: 'bool' = True) -> 'dict[str, Any]'`。
- **输入：** 必填：`merged_calls`, `output_directory`；其余参数均为仅关键字参数且可选（`None` 或 `False` 表示使用工具自己的默认值）。
- **输出：** 含 `status, backend, operation, version, command, output_directory, outputs, log_tail, provenance` 的字典，另有该函数特有的键（见英文版 Output）。

## 7) `dnakit.alignment.minigraph_p11_cli_functions.convert_to_stable_gaf()`

- **作用：** 把图比对（GAF）转换为线性序列的稳定坐标（mgutils.js stableGaf）。
- **处理方法：** 运行固定的参数列表（不经过 shell）；参数含义、取值范围和工具默认值见英文版同名小节的 Method 说明。
- **API：** `dnakit.alignment.minigraph_p11_cli_functions.convert_to_stable_gaf(graph: 'str | Path', alignments: 'str | Path', output_directory: 'str | Path', *, output_name: 'str' = 'stable.gaf', mgutils_path: 'str | None' = None, k8_path: 'str | None' = None, timeout_seconds: 'float' = 3600.0, verbose: 'bool' = True) -> 'dict[str, Any]'`。
- **输入：** 必填：`graph`, `alignments`, `output_directory`；其余参数均为仅关键字参数且可选（`None` 或 `False` 表示使用工具自己的默认值）。
- **输出：** 含 `status, backend, operation, version, command, output_directory, outputs, log_tail, provenance` 的字典，另有该函数特有的键（见英文版 Output）。

## 8) `dnakit.alignment.minigraph_p11_cli_functions.count_gene_copies()`

- **作用：** 由基因序列的图比对统计组装中的基因拷贝数（mgutils.js genecopy）。
- **处理方法：** 运行固定的参数列表（不经过 shell）；参数含义、取值范围和工具默认值见英文版同名小节的 Method 说明。
- **API：** `dnakit.alignment.minigraph_p11_cli_functions.count_gene_copies(alignments: 'str | Path', gene_bed: 'str | Path', output_directory: 'str | Path', *, output_name: 'str' = 'gene_copy.tsv', min_coverage: 'float | None' = None, min_relative_coverage: 'float | None' = None, mgutils_path: 'str | None' = None, k8_path: 'str | None' = None, timeout_seconds: 'float' = 3600.0, verbose: 'bool' = True) -> 'dict[str, Any]'`。
- **输入：** 必填：`alignments`, `gene_bed`, `output_directory`；其余参数均为仅关键字参数且可选（`None` 或 `False` 表示使用工具自己的默认值）。
- **输出：** 含 `status, backend, operation, version, command, output_directory, outputs, log_tail, provenance` 的字典，另有该函数特有的键（见英文版 Output）。

## 9) `dnakit.alignment.minigraph_p11_cli_functions.evaluate_sv_calls()`

- **作用：** 用真集评估结构变异调用，或比较多个调用集（gafcall.js eval）。
- **处理方法：** 运行固定的参数列表（不经过 shell）；参数含义、取值范围和工具默认值见英文版同名小节的 Method 说明。
- **API：** `dnakit.alignment.minigraph_p11_cli_functions.evaluate_sv_calls(vcf_files: 'Sequence[str | Path]', output_directory: 'str | Path', *, output_name: 'str' = 'eval.txt', confident_regions: 'str | Path | None' = None, min_sv_length: 'int | None' = None, window_size: 'int | None' = None, read_length_ratio: 'float | None' = None, min_length_ratio: 'float | None' = None, min_allele_frequency: 'float | None' = None, ignore_filter: 'bool' = False, check_genotype: 'bool' = False, print_errors: 'bool' = False, gafcall_path: 'str | None' = None, k8_path: 'str | None' = None, timeout_seconds: 'float' = 3600.0, verbose: 'bool' = True) -> 'dict[str, Any]'`。
- **输入：** 必填：`vcf_files`, `output_directory`；其余参数均为仅关键字参数且可选（`None` 或 `False` 表示使用工具自己的默认值）。
- **输出：** 含 `status, backend, operation, version, command, output_directory, outputs, log_tail, provenance` 的字典，另有该函数特有的键（见英文版 Output）。

## 10) `dnakit.alignment.minigraph_p11_cli_functions.extract_gaf_segment_path()`

- **作用：** 从 GAF 比对中提取两个图片段之间的 contig 区间（mgutils.js extractseg）。
- **处理方法：** 运行固定的参数列表（不经过 shell）；参数含义、取值范围和工具默认值见英文版同名小节的 Method 说明。
- **API：** `dnakit.alignment.minigraph_p11_cli_functions.extract_gaf_segment_path(first_segment: 'str', second_segment: 'str', alignments: 'Sequence[str | Path]', output_directory: 'str | Path', *, output_name: 'str' = 'segment_path.bed', min_alignment_length: 'int' = 100000, from_end: 'bool' = False, mgutils_path: 'str | None' = None, k8_path: 'str | None' = None, timeout_seconds: 'float' = 3600.0, verbose: 'bool' = True) -> 'dict[str, Any]'`。
- **输入：** 必填：`first_segment`, `second_segment`, `alignments`, `output_directory`；其余参数均为仅关键字参数且可选（`None` 或 `False` 表示使用工具自己的默认值）。
- **输出：** 含 `status, backend, operation, version, command, output_directory, outputs, log_tail, provenance` 的字典，另有该函数特有的键（见英文版 Output）。

## 11) `dnakit.alignment.minigraph_p11_cli_functions.extract_long_indels()`

- **作用：** 从 reads 的稳定坐标 GAF 比对中提取长 INDEL 和断点（gafcall.js extract）。
- **处理方法：** 运行固定的参数列表（不经过 shell）；参数含义、取值范围和工具默认值见英文版同名小节的 Method 说明。
- **API：** `dnakit.alignment.minigraph_p11_cli_functions.extract_long_indels(stable_gaf: 'str | Path', output_directory: 'str | Path', *, output_name: 'str' = 'extract.tsv', sample_name: 'str | None' = None, min_indel_length: 'int | None' = None, max_indels_per_10kb: 'int | None' = None, min_mapq: 'int | None' = None, min_mapq_ends: 'int | None' = None, min_alignment_length_ends: 'int | None' = None, min_alignment_length_middle: 'int | None' = None, polya_penalty: 'int | None' = None, centromere_bed: 'str | Path | None' = None, gafcall_path: 'str | None' = None, k8_path: 'str | None' = None, timeout_seconds: 'float' = 3600.0, verbose: 'bool' = True) -> 'dict[str, Any]'`。
- **输入：** 必填：`stable_gaf`, `output_directory`；其余参数均为仅关键字参数且可选（`None` 或 `False` 表示使用工具自己的默认值）。
- **输出：** 含 `status, backend, operation, version, command, output_directory, outputs, log_tail, provenance` 的字典，另有该函数特有的键（见英文版 Output）。

## 12) `dnakit.alignment.minigraph_p11_cli_functions.filter_sv_extract_by_type()`

- **作用：** 保留一个 SV 提取文件中其 read 在另一文件里有相同 SV 类型的记录（gafcall.js join）。
- **处理方法：** 运行固定的参数列表（不经过 shell）；参数含义、取值范围和工具默认值见英文版同名小节的 Method 说明。
- **API：** `dnakit.alignment.minigraph_p11_cli_functions.filter_sv_extract_by_type(filter_extract: 'str | Path', extract: 'str | Path', output_directory: 'str | Path', *, output_name: 'str' = 'joined.tsv', gafcall_path: 'str | None' = None, k8_path: 'str | None' = None, timeout_seconds: 'float' = 3600.0, verbose: 'bool' = True) -> 'dict[str, Any]'`。
- **输入：** 必填：`filter_extract`, `extract`, `output_directory`；其余参数均为仅关键字参数且可选（`None` 或 `False` 表示使用工具自己的默认值）。
- **输出：** 含 `status, backend, operation, version, command, output_directory, outputs, log_tail, provenance` 的字典，另有该函数特有的键（见英文版 Output）。

## 13) `dnakit.alignment.minigraph_p11_cli_functions.graph_segment_coverage()`

- **作用：** 把 reads 重新映射到图上并报告每个图片段的覆盖度（minigraph --cov）。
- **处理方法：** 运行固定的参数列表（不经过 shell）；参数含义、取值范围和工具默认值见英文版同名小节的 Method 说明。
- **API：** `dnakit.alignment.minigraph_p11_cli_functions.graph_segment_coverage(graph: 'str | Path', reads: 'Sequence[str | Path]', output_directory: 'str | Path', *, output_name: 'str' = 'coverage.gfa', preset: 'str | None' = None, base_alignment: 'bool' = True, min_aligned_length: 'int | None' = None, min_mapq: 'int | None' = None, kmer_size: 'int | None' = None, window_size: 'int | None' = None, threads: 'int' = 4, minigraph_path: 'str | None' = None, timeout_seconds: 'float' = 3600.0, verbose: 'bool' = True) -> 'dict[str, Any]'`。
- **输入：** 必填：`graph`, `reads`, `output_directory`；其余参数均为仅关键字参数且可选（`None` 或 `False` 表示使用工具自己的默认值）。
- **输出：** 含 `status, backend, operation, version, command, output_directory, outputs, log_tail, provenance` 的字典，另有该函数特有的键（见英文版 Output）。

## 14) `dnakit.alignment.minigraph_p11_cli_functions.map_sequences_to_graph()`

- **作用：** 把 DNA 序列映射到序列图或线性参考（minigraph）。
- **处理方法：** 运行固定的参数列表（不经过 shell）；参数含义、取值范围和工具默认值见英文版同名小节的 Method 说明。
- **API：** `dnakit.alignment.minigraph_p11_cli_functions.map_sequences_to_graph(graph: 'str | Path', queries: 'Sequence[str | Path]', output_directory: 'str | Path', *, output_name: 'str' = 'mapping.gaf', preset: 'str | None' = None, base_alignment: 'bool' = False, kmer_size: 'int | None' = None, window_size: 'int | None' = None, repetitive_fraction: 'float | None' = None, occurrence_range: 'int | Sequence[int] | None' = None, divergence: 'float | None' = None, chain_gap: 'int | None' = None, bandwidth: 'int | Sequence[int] | None' = None, min_minimizers: 'int | Sequence[int] | None' = None, min_chain_score: 'int | Sequence[int] | None' = None, min_secondary_ratio: 'float | None' = None, max_secondary: 'int | None' = None, skip_self_diagonal: 'bool' = False, linear_chains: 'bool' = False, vertex_coordinates: 'bool' = False, threads: 'int' = 4, minigraph_path: 'str | None' = None, timeout_seconds: 'float' = 3600.0, verbose: 'bool' = True) -> 'dict[str, Any]'`。
- **输入：** 必填：`graph`, `queries`, `output_directory`；其余参数均为仅关键字参数且可选（`None` 或 `False` 表示使用工具自己的默认值）。
- **输出：** 含 `status, backend, operation, version, command, output_directory, outputs, log_tail, provenance` 的字典，另有该函数特有的键（见英文版 Output）。

## 15) `dnakit.alignment.minigraph_p11_cli_functions.merge_bubble_calls()`

- **作用：** 把各样本的 --call BED 合并为多样本表（mgutils.js merge）。
- **处理方法：** 运行固定的参数列表（不经过 shell）；参数含义、取值范围和工具默认值见英文版同名小节的 Method 说明。
- **API：** `dnakit.alignment.minigraph_p11_cli_functions.merge_bubble_calls(call_files: 'Sequence[str | Path]', output_directory: 'str | Path', *, output_name: 'str' = 'merged.txt', samples: 'Sequence[str] | None' = None, annotation: 'str | Path | None' = None, mgutils_path: 'str | None' = None, k8_path: 'str | None' = None, timeout_seconds: 'float' = 3600.0, verbose: 'bool' = True) -> 'dict[str, Any]'`。
- **输入：** 必填：`call_files`, `output_directory`；其余参数均为仅关键字参数且可选（`None` 或 `False` 表示使用工具自己的默认值）。
- **输出：** 含 `status, backend, operation, version, command, output_directory, outputs, log_tail, provenance` 的字典，另有该函数特有的键（见英文版 Output）。

## 16) `dnakit.alignment.minigraph_p11_cli_functions.merge_sv_evidence()`

- **作用：** 把 read 级 INDEL/断点证据合并为结构变异调用（gafcall.js merge）。
- **处理方法：** 运行固定的参数列表（不经过 shell）；参数含义、取值范围和工具默认值见英文版同名小节的 Method 说明。
- **API：** `dnakit.alignment.minigraph_p11_cli_functions.merge_sv_evidence(extract_files: 'Sequence[str | Path]', output_directory: 'str | Path', *, output_name: 'str' = 'calls.tsv', min_read_count: 'int | None' = None, min_read_count_per_strand: 'int | None' = None, window_size: 'int | None' = None, max_length_difference: 'float | None' = None, min_centromere_distance: 'int | None' = None, min_retrotransposon_length: 'int | None' = None, min_retrotransposon_reads: 'int | None' = None, max_alleles: 'int | None' = None, max_reads_compared: 'int | None' = None, gafcall_path: 'str | None' = None, k8_path: 'str | None' = None, timeout_seconds: 'float' = 3600.0, verbose: 'bool' = True) -> 'dict[str, Any]'`。
- **输入：** 必填：`extract_files`, `output_directory`；其余参数均为仅关键字参数且可选（`None` 或 `False` 表示使用工具自己的默认值）。
- **输出：** 含 `status, backend, operation, version, command, output_directory, outputs, log_tail, provenance` 的字典，另有该函数特有的键（见英文版 Output）。

## 17) `dnakit.alignment.minigraph_p11_cli_functions.rename_fasta_sequences()`

- **作用：** 给 FASTA 序列名加前缀（PanSN 约定；mgutils.js renamefa）。
- **处理方法：** 运行固定的参数列表（不经过 shell）；参数含义、取值范围和工具默认值见英文版同名小节的 Method 说明。
- **API：** `dnakit.alignment.minigraph_p11_cli_functions.rename_fasta_sequences(fasta: 'str | Path', prefix: 'str', output_directory: 'str | Path', *, output_name: 'str' = 'renamed.fa', delimiter: 'str' = '#', mgutils_path: 'str | None' = None, k8_path: 'str | None' = None, timeout_seconds: 'float' = 3600.0, verbose: 'bool' = True) -> 'dict[str, Any]'`。
- **输入：** 必填：`fasta`, `prefix`, `output_directory`；其余参数均为仅关键字参数且可选（`None` 或 `False` 表示使用工具自己的默认值）。
- **输出：** 含 `status, backend, operation, version, command, output_directory, outputs, log_tail, provenance` 的字典，另有该函数特有的键（见英文版 Output）。

## 18) `dnakit.alignment.minigraph_p11_cli_functions.segment_allele_frequency()`

- **作用：** 由合并的气泡结果计算每个图片段的等位频率（mgutils.js segfreq）。
- **处理方法：** 运行固定的参数列表（不经过 shell）；参数含义、取值范围和工具默认值见英文版同名小节的 Method 说明。
- **API：** `dnakit.alignment.minigraph_p11_cli_functions.segment_allele_frequency(segment_bed: 'str | Path', merged_calls: 'str | Path', output_directory: 'str | Path', *, output_name: 'str' = 'segfreq.tsv', bubble_bed: 'str | Path | None' = None, min_frequency: 'float | None' = None, mgutils_path: 'str | None' = None, k8_path: 'str | None' = None, timeout_seconds: 'float' = 3600.0, verbose: 'bool' = True) -> 'dict[str, Any]'`。
- **输入：** 必填：`segment_bed`, `merged_calls`, `output_directory`；其余参数均为仅关键字参数且可选（`None` 或 `False` 表示使用工具自己的默认值）。
- **输出：** 含 `status, backend, operation, version, command, output_directory, outputs, log_tail, provenance` 的字典，另有该函数特有的键（见英文版 Output）。

## 19) `dnakit.alignment.minigraph_p11_cli_functions.summarize_insert_annotation()`

- **作用：** 把插入序列注释汇总为重复类别表（mgutils.js anno2tbl）。
- **处理方法：** 运行固定的参数列表（不经过 shell）；参数含义、取值范围和工具默认值见英文版同名小节的 Method 说明。
- **API：** `dnakit.alignment.minigraph_p11_cli_functions.summarize_insert_annotation(annotation: 'str | Path', output_directory: 'str | Path', *, output_name: 'str' = 'anno_table.tsv', mgutils_path: 'str | None' = None, k8_path: 'str | None' = None, timeout_seconds: 'float' = 3600.0, verbose: 'bool' = True) -> 'dict[str, Any]'`。
- **输入：** 必填：`annotation`, `output_directory`；其余参数均为仅关键字参数且可选（`None` 或 `False` 表示使用工具自己的默认值）。
- **输出：** 含 `status, backend, operation, version, command, output_directory, outputs, log_tail, provenance` 的字典，另有该函数特有的键（见英文版 Output）。

## 20) `dnakit.alignment.minigraph_p11_cli_functions.view_sv_calls()`

- **作用：** 把 VCF 中的结构变异重新格式化为带长 INDEL 计数的表（gafcall.js view）。
- **处理方法：** 运行固定的参数列表（不经过 shell）；参数含义、取值范围和工具默认值见英文版同名小节的 Method 说明。
- **API：** `dnakit.alignment.minigraph_p11_cli_functions.view_sv_calls(vcf: 'str | Path', output_directory: 'str | Path', *, output_name: 'str' = 'sv.tsv', min_length: 'int | None' = None, regions: 'str | Path | None' = None, ignore_filter: 'bool' = False, check_genotype: 'bool' = False, count_long: 'bool' = False, gafcall_path: 'str | None' = None, k8_path: 'str | None' = None, timeout_seconds: 'float' = 3600.0, verbose: 'bool' = True) -> 'dict[str, Any]'`。
- **输入：** 必填：`vcf`, `output_directory`；其余参数均为仅关键字参数且可选（`None` 或 `False` 表示使用工具自己的默认值）。
- **输出：** 含 `status, backend, operation, version, command, output_directory, outputs, log_tail, provenance` 的字典，另有该函数特有的键（见英文版 Output）。

