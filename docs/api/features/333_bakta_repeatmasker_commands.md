# Bakta 与 RepeatMasker：基因组注释、重复掩蔽与重复注释工具

`dnakit.annotation` 中的这一组函数是对 [Bakta](https://github.com/oschwengers/bakta) 和 [RepeatMasker](https://www.repeatmasker.org/) 的命令行封装。Bakta：`annotate_bakta_features` 注释细菌基因组并可选择检测步骤（CDS 及其数据库功能注释、AMR 基因与专家蛋白、假基因、sORF、tRNA、tmRNA、rRNA、ncRNA、ncRNA 区域、CRISPR、组装缺口、oriC/oriV/oriT；跳过的步骤对应 `--skip-*`），`export_bakta_annotation` 把 Bakta JSON 重新写成全部格式（TSV、推断 TSV、hypothetical、GenBank、EMBL、GFF3、核酸/蛋白 FASTA、汇总），`plot_bakta_genome` 绘制环形图（含 COG 图）。RepeatMasker：`mask_repeats`、`clip_bacterial_insertion_sequences`、`process_repeat_hits`、`date_repeats`、`annotate_segmental_duplications`、`mask_repeat_proteins`、`mask_tandem_repeats` 运行相应程序，`util/` 脚本对应 `convert_repeats_to_gff3`、`convert_repeats_to_bed`、`export_repeats_ucsc_tables`、`build_repeat_track_hub`、`extract_repeat_sequences`、`mask_sequences_from_annotation`、`summarize_repeat_annotation`、`recalculate_repeat_divergence`、`plot_repeat_landscape`、`plot_duplicon_map`、`combine_repeat_annotations`。未提供：`wublastToCrossmatch.pl`（重排已停用的 WU-BLAST 引擎输出）和 `addRepBase.pl`（需要有许可的 RepBase RepeatMasker Edition）。

测试版本为 Bakta 1.12.2（GPL-3.0）和 RepeatMasker 4.2.6（Open Software License），均作为独立程序启动。用 `conda install -c bioconda bakta` 安装 Bakta（并 `bakta_db download` 下载数据库），用 `conda install -c bioconda repeatmasker` 安装 RepeatMasker（然后 `perl configure`），把程序放到 PATH 上或传入 `bakta_path`、`bakta_io_path`、`bakta_plot_path`、`repeatmasker_dir`。每个函数都写入新的 `output_directory`；RepeatMasker 程序把结果写在输入旁边，因此输入会被复制到其中。完整的 `bakta` 注释需要数据库和多个外部工具，仅验证到参数解析；`bakta_io`、`bakta_plot` 和 RepeatMasker 函数（需要 FamDB 的 `clip_bacterial_insertion_sequences`、`date_repeats`、`annotate_segmental_duplications`、`mask_repeat_proteins` 除外）均真实运行过。示例使用 `verbose=False`。

## 示例数据与示例

**示例文件：** `example.out`

```text
   SW   perc perc perc  query     position in query    matching repeat       position in repeat
score   div. del. ins.  sequence  begin end   (left)   repeat   class/family begin  end    (left)  ID

  734   27.6  9.9  4.0  seq1       1483  1806  (326) + L1M5     LINE/L1      5368   5710  (414)   1
 2134   10.7  3.3  1.3  seq1       1809  2111   (21) + AluSx    SINE/Alu        3    311    (1)   2
```

**示例代码：**

```python
from pathlib import Path
import tempfile
from dnakit.annotation import repeatmasker_p11_cli_functions as rm

work = Path(tempfile.mkdtemp())
gff = rm.convert_repeats_to_gff3("example.out", work / "gff", verbose=False)
print(gff["features"])
print(Path(gff["outputs"]["repeats.gff3"]).read_text().splitlines()[2])
bed = rm.convert_repeats_to_bed("example.out", work / "bed", verbose=False)
print(bed["intervals"])
```

**示例输出：**

```text
2
seq1	RepeatMasker	dispersed_repeat	1483	1806	734	+	.	ID=1;Target=L1M5 5368 5710
2
```

**示例代码：**

```python
from pathlib import Path
import json, random, tempfile
from dnakit.annotation import bakta_p11_cli_functions as bk

rng = random.Random(2)
seq = "".join(rng.choice("ACGT") for _ in range(3000))
cds = {"sequence": "c1", "strand": "+", "nt": seq[100:166], "aa": "MKTAYIAKQRQISFVKSHFSRQ", "aa_digest": None,
       "start_type": "ATG", "rbs_motif": None, "frame": 2, "functions": [], "genes": []}
data = {
    "genome": {"genus": "E", "species": "coli", "strain": "", "complete": False, "gram": "?", "translation_table": 11},
    "stats": {"no_sequences": 1, "size": 3000, "gc": 0.5, "n_ratio": 0.0, "n50": 3000, "n90": 3000, "coding_ratio": 0.1},
    "features": [
        {**cds, "type": "cds", "start": 101, "stop": 166, "gene": "abcD", "product": "test protein", "locus": "T_00005", "db_xrefs": ["COG:J"]},
        {"type": "tRNA", "sequence": "c1", "start": 501, "stop": 576, "strand": "+", "gene": "tRNA-Ala", "product": "tRNA-Ala(tgc)",
         "locus": "T_00015", "nt": seq[500:576], "amino_acid": "Ala", "anti_codon": "tgc", "score": 80.5, "db_xrefs": ["SO:0001272"]},
    ],
    "sequences": [{"id": "c1", "description": "", "nt": seq, "length": 3000, "complete": False, "type": "contig",
                   "topology": "linear", "name": None, "orig_id": "c1", "orig_description": ""}],
    "version": {"bakta": "1.12.2", "db": {"type": "light", "version": "6.0", "doi": "10.5281/zenodo.0"}},
    "run": {"start": "2026-10-20 10:00:00", "end": "2026-10-20 10:01:00"},
}
work = Path(tempfile.mkdtemp())
(work / "demo.json").write_text(json.dumps(data))
export = bk.export_bakta_annotation(work / "demo.json", work / "export", verbose=False)
print(export["version"], sorted(k for k in export["outputs"] if k in ("gbff", "embl", "gff3", "tsv", "faa", "ffn", "txt")))
print(export["summary"]["CDSs"], export["summary"]["tRNAs"])
plot = bk.plot_bakta_genome(work / "demo.json", work / "plot", plot_type="cog", size=4, dpi=150, verbose=False)
print(sorted(k for k in plot["outputs"] if k in ("png", "svg")))
```

**示例输出：**

```text
1.12.2 ['embl', 'faa', 'ffn', 'gbff', 'gff3', 'tsv', 'txt']
1 1
['png', 'svg']
```


## 1) `dnakit.annotation.bakta_p11_cli_functions.annotate_bakta_features()`

- **作用：** 用 Bakta 注释细菌基因组、MAG 或质粒，并可选择要检测的特征类型。
- **处理方法：** 运行固定的参数列表（不经过 shell）；参数含义、取值范围和工具默认值见英文版同名小节的 Method 说明。
- **API：** `dnakit.annotation.bakta_p11_cli_functions.annotate_bakta_features(genome: 'str | Path', database: 'str | Path', output_directory: 'str | Path', *, features: 'Sequence[str] | None' = None, prefix: 'str' = 'annotation', min_contig_length: 'int | None' = None, genus: 'str | None' = None, species: 'str | None' = None, strain: 'str | None' = None, plasmid: 'str | None' = None, complete: 'bool' = False, prodigal_training_file: 'str | Path | None' = None, translation_table: 'int | None' = None, gram: 'str | None' = None, locus: 'str | None' = None, locus_tag: 'str | None' = None, locus_tag_increment: 'int | None' = None, keep_contig_headers: 'bool' = False, compliant: 'bool' = False, replicons: 'str | Path | None' = None, regions: 'str | Path | None' = None, trusted_proteins: 'str | Path | None' = None, trusted_hmms: 'str | Path | None' = None, metagenome: 'bool' = False, partial: 'bool' = False, skip_overlap_filters: 'bool' = False, plot: 'bool' = False, threads: 'int | None' = None, bakta_path: 'str | None' = None, timeout_seconds: 'float' = 86400.0, verbose: 'bool' = True) -> 'dict[str, Any]'`。
- **输入：** 必填：`genome`, `database`, `output_directory`；其余参数均为仅关键字参数且可选（`None` 或 `False` 表示使用工具自己的默认值）。
- **输出：** 含 `status, backend, operation, version, command, output_directory, outputs, log_tail, provenance` 的字典，另有该函数特有的键（见英文版 Output）。

## 2) `dnakit.annotation.bakta_p11_cli_functions.export_bakta_annotation()`

- **作用：** 把 Bakta 的 JSON 注释重新写成所有输出格式，无需重新运行注释。
- **处理方法：** 运行固定的参数列表（不经过 shell）；参数含义、取值范围和工具默认值见英文版同名小节的 Method 说明。
- **API：** `dnakit.annotation.bakta_p11_cli_functions.export_bakta_annotation(annotation_json: 'str | Path', output_directory: 'str | Path', *, prefix: 'str | None' = None, bakta_io_path: 'str | None' = None, timeout_seconds: 'float' = 3600.0, verbose: 'bool' = True) -> 'dict[str, Any]'`。
- **输入：** 必填：`annotation_json`, `output_directory`；其余参数均为仅关键字参数且可选（`None` 或 `False` 表示使用工具自己的默认值）。
- **输出：** 含 `status, backend, operation, version, command, output_directory, outputs, log_tail, provenance` 的字典，另有该函数特有的键（见英文版 Output）。

## 3) `dnakit.annotation.bakta_p11_cli_functions.plot_bakta_genome()`

- **作用：** 由 Bakta JSON 注释绘制环形基因组图（特征类型或 COG 功能类别）。
- **处理方法：** 运行固定的参数列表（不经过 shell）；参数含义、取值范围和工具默认值见英文版同名小节的 Method 说明。
- **API：** `dnakit.annotation.bakta_p11_cli_functions.plot_bakta_genome(annotation_json: 'str | Path', output_directory: 'str | Path', *, plot_type: 'str' = 'features', prefix: 'str | None' = None, sequences: 'Sequence[str] | None' = None, config: 'str | Path | None' = None, label: 'str | None' = None, size: 'int | None' = None, dpi: 'int | None' = None, bakta_plot_path: 'str | None' = None, timeout_seconds: 'float' = 3600.0, verbose: 'bool' = True) -> 'dict[str, Any]'`。
- **输入：** 必填：`annotation_json`, `output_directory`；其余参数均为仅关键字参数且可选（`None` 或 `False` 表示使用工具自己的默认值）。
- **输出：** 含 `status, backend, operation, version, command, output_directory, outputs, log_tail, provenance` 的字典，另有该函数特有的键（见英文版 Output）。

## 4) `dnakit.annotation.repeatmasker_p11_cli_functions.annotate_segmental_duplications()`

- **作用：** 用重复片段库搜索并注释片段重复（DupMasker）。
- **处理方法：** 运行固定的参数列表（不经过 shell）；参数含义、取值范围和工具默认值见英文版同名小节的 Method 说明。
- **API：** `dnakit.annotation.repeatmasker_p11_cli_functions.annotate_segmental_duplications(fasta: 'str | Path', output_directory: 'str | Path', *, repeatmasker_out: 'str | Path | None' = None, max_divergence: 'float | None' = None, max_gap_width: 'int | None' = None, force_search: 'bool' = False, alignments: 'bool' = False, gff_output: 'bool' = False, search_engine: 'str | None' = None, parallel: 'int | None' = None, seeds_only: 'bool' = False, repeatmasker_dir: 'str | Path | None' = None, timeout_seconds: 'float' = 86400.0, verbose: 'bool' = True) -> 'dict[str, Any]'`。
- **输入：** 必填：`fasta`, `output_directory`；其余参数均为仅关键字参数且可选（`None` 或 `False` 表示使用工具自己的默认值）。
- **输出：** 含 `status, backend, operation, version, command, output_directory, outputs, log_tail, provenance` 的字典，另有该函数特有的键（见英文版 Output）。

## 5) `dnakit.annotation.repeatmasker_p11_cli_functions.build_repeat_track_hub()`

- **作用：** 由重复注释构建 UCSC bigRmsk track hub（rmToTrackHub.pl）。
- **处理方法：** 运行固定的参数列表（不经过 shell）；参数含义、取值范围和工具默认值见英文版同名小节的 Method 说明。
- **API：** `dnakit.annotation.repeatmasker_p11_cli_functions.build_repeat_track_hub(repeat_annotation: 'str | Path', output_directory: 'str | Path', *, genome: 'str', hub_name: 'str', alignments: 'str | Path | None' = None, chrom_sizes: 'str | Path | None' = None, only_bigbed: 'bool' = False, repeatmasker_dir: 'str | Path | None' = None, timeout_seconds: 'float' = 3600.0, verbose: 'bool' = True) -> 'dict[str, Any]'`。
- **输入：** 必填：`repeat_annotation`, `output_directory`；其余参数均为仅关键字参数且可选（`None` 或 `False` 表示使用工具自己的默认值）。
- **输出：** 含 `status, backend, operation, version, command, output_directory, outputs, log_tail, provenance` 的字典，另有该函数特有的键（见英文版 Output）。

## 6) `dnakit.annotation.repeatmasker_p11_cli_functions.clip_bacterial_insertion_sequences()`

- **作用：** 剪除 DNA 中的细菌插入序列（IS）元件（RepeatMasker -is_only）。
- **处理方法：** 运行固定的参数列表（不经过 shell）；参数含义、取值范围和工具默认值见英文版同名小节的 Method 说明。
- **API：** `dnakit.annotation.repeatmasker_p11_cli_functions.clip_bacterial_insertion_sequences(fasta: 'str | Path', output_directory: 'str | Path', *, repeatmasker_dir: 'str | Path | None' = None, timeout_seconds: 'float' = 86400.0, verbose: 'bool' = True) -> 'dict[str, Any]'`。
- **输入：** 必填：`fasta`, `output_directory`；其余参数均为仅关键字参数且可选（`None` 或 `False` 表示使用工具自己的默认值）。
- **输出：** 含 `status, backend, operation, version, command, output_directory, outputs, log_tail, provenance` 的字典，另有该函数特有的键（见英文版 Output）。

## 7) `dnakit.annotation.repeatmasker_p11_cli_functions.combine_repeat_annotations()`

- **作用：** 把两套 RepeatMasker 注释（.out/.align）合并为一套排序后的结果（combineRMFiles.pl）。
- **处理方法：** 运行固定的参数列表（不经过 shell）；参数含义、取值范围和工具默认值见英文版同名小节的 Method 说明。
- **API：** `dnakit.annotation.repeatmasker_p11_cli_functions.combine_repeat_annotations(first_annotation: 'str | Path', first_alignments: 'str | Path', second_annotation: 'str | Path', second_alignments: 'str | Path', output_directory: 'str | Path', *, output_prefix: 'str' = 'combined', repeatmasker_dir: 'str | Path | None' = None, timeout_seconds: 'float' = 3600.0, verbose: 'bool' = True) -> 'dict[str, Any]'`。
- **输入：** 必填：`first_annotation`, `first_alignments`, `second_annotation`, `second_alignments`, `output_directory`；其余参数均为仅关键字参数且可选（`None` 或 `False` 表示使用工具自己的默认值）。
- **输出：** 含 `status, backend, operation, version, command, output_directory, outputs, log_tail, provenance` 的字典，另有该函数特有的键（见英文版 Output）。

## 8) `dnakit.annotation.repeatmasker_p11_cli_functions.convert_repeats_to_bed()`

- **作用：** 把重复注释转换为 BED，可按类别、长度和分歧度过滤（RM2Bed.py）。
- **处理方法：** 运行固定的参数列表（不经过 shell）；参数含义、取值范围和工具默认值见英文版同名小节的 Method 说明。
- **API：** `dnakit.annotation.repeatmasker_p11_cli_functions.convert_repeats_to_bed(repeat_file: 'str | Path', output_directory: 'str | Path', *, output_prefix: 'str | None' = None, split_by: 'str | None' = None, sort_by: 'str | None' = None, min_length: 'int | None' = None, min_hits: 'int | None' = None, max_divergence: 'float | None' = None, min_divergence: 'float | None' = None, overlap_resolution: 'str | None' = None, repeatmasker_dir: 'str | Path | None' = None, timeout_seconds: 'float' = 3600.0, verbose: 'bool' = True) -> 'dict[str, Any]'`。
- **输入：** 必填：`repeat_file`, `output_directory`；其余参数均为仅关键字参数且可选（`None` 或 `False` 表示使用工具自己的默认值）。
- **输出：** 含 `status, backend, operation, version, command, output_directory, outputs, log_tail, provenance` 的字典，另有该函数特有的键（见英文版 Output）。

## 9) `dnakit.annotation.repeatmasker_p11_cli_functions.convert_repeats_to_gff3()`

- **作用：** 把重复注释（.out）转换为 GFF3（rmOutToGFF3.pl）。
- **处理方法：** 运行固定的参数列表（不经过 shell）；参数含义、取值范围和工具默认值见英文版同名小节的 Method 说明。
- **API：** `dnakit.annotation.repeatmasker_p11_cli_functions.convert_repeats_to_gff3(repeat_annotation: 'str | Path', output_directory: 'str | Path', *, output_name: 'str' = 'repeats.gff3', repeatmasker_dir: 'str | Path | None' = None, timeout_seconds: 'float' = 3600.0, verbose: 'bool' = True) -> 'dict[str, Any]'`。
- **输入：** 必填：`repeat_annotation`, `output_directory`；其余参数均为仅关键字参数且可选（`None` 或 `False` 表示使用工具自己的默认值）。
- **输出：** 含 `status, backend, operation, version, command, output_directory, outputs, log_tail, provenance` 的字典，另有该函数特有的键（见英文版 Output）。

## 10) `dnakit.annotation.repeatmasker_p11_cli_functions.date_repeats()`

- **作用：** 标注重复是谱系特异还是祖先重复，并可掩蔽谱系特异重复（DateRepeats）。
- **处理方法：** 运行固定的参数列表（不经过 shell）；参数含义、取值范围和工具默认值见英文版同名小节的 Method 说明。
- **API：** `dnakit.annotation.repeatmasker_p11_cli_functions.date_repeats(repeat_annotation: 'str | Path', output_directory: 'str | Path', *, query_species: 'str', comparison_species: 'Sequence[str]', sequence_file: 'str | Path | None' = None, mask_species: 'str | None' = None, mask_uncertain: 'bool' = False, mask_low_complexity: 'bool' = True, repeatmasker_dir: 'str | Path | None' = None, timeout_seconds: 'float' = 86400.0, verbose: 'bool' = True) -> 'dict[str, Any]'`。
- **输入：** 必填：`repeat_annotation`, `output_directory`；其余参数均为仅关键字参数且可选（`None` 或 `False` 表示使用工具自己的默认值）。
- **输出：** 含 `status, backend, operation, version, command, output_directory, outputs, log_tail, provenance` 的字典，另有该函数特有的键（见英文版 Output）。

## 11) `dnakit.annotation.repeatmasker_p11_cli_functions.export_repeats_ucsc_tables()`

- **作用：** 把重复注释导出为 UCSC rmsk 表（rmToUCSCTables.pl）。
- **处理方法：** 运行固定的参数列表（不经过 shell）；参数含义、取值范围和工具默认值见英文版同名小节的 Method 说明。
- **API：** `dnakit.annotation.repeatmasker_p11_cli_functions.export_repeats_ucsc_tables(repeat_annotation: 'str | Path', output_directory: 'str | Path', *, alignments: 'str | Path | None' = None, repeatmasker_dir: 'str | Path | None' = None, timeout_seconds: 'float' = 3600.0, verbose: 'bool' = True) -> 'dict[str, Any]'`。
- **输入：** 必填：`repeat_annotation`, `output_directory`；其余参数均为仅关键字参数且可选（`None` 或 `False` 表示使用工具自己的默认值）。
- **输出：** 含 `status, backend, operation, version, command, output_directory, outputs, log_tail, provenance` 的字典，另有该函数特有的键（见英文版 Output）。

## 12) `dnakit.annotation.repeatmasker_p11_cli_functions.extract_repeat_sequences()`

- **作用：** 把每个注释的重复元件序列提取为 FASTA（rmOut2Fasta.pl）。
- **处理方法：** 运行固定的参数列表（不经过 shell）；参数含义、取值范围和工具默认值见英文版同名小节的 Method 说明。
- **API：** `dnakit.annotation.repeatmasker_p11_cli_functions.extract_repeat_sequences(fasta: 'str | Path', repeat_annotation: 'str | Path', output_directory: 'str | Path', *, output_name: 'str' = 'repeats.fa', repeatmasker_dir: 'str | Path | None' = None, timeout_seconds: 'float' = 3600.0, verbose: 'bool' = True) -> 'dict[str, Any]'`。
- **输入：** 必填：`fasta`, `repeat_annotation`, `output_directory`；其余参数均为仅关键字参数且可选（`None` 或 `False` 表示使用工具自己的默认值）。
- **输出：** 含 `status, backend, operation, version, command, output_directory, outputs, log_tail, provenance` 的字典，另有该函数特有的键（见英文版 Output）。

## 13) `dnakit.annotation.repeatmasker_p11_cli_functions.mask_repeat_proteins()`

- **作用：** 通过翻译搜索掩蔽 DNA 中重复来源的蛋白编码片段（RepeatProteinMask）。
- **处理方法：** 运行固定的参数列表（不经过 shell）；参数含义、取值范围和工具默认值见英文版同名小节的 Method 说明。
- **API：** `dnakit.annotation.repeatmasker_p11_cli_functions.mask_repeat_proteins(fasta: 'str | Path', output_directory: 'str | Path', *, search_engine: 'str | None' = None, p_value: 'float | None' = None, min_score: 'int | None' = None, word_size: 'int | None' = None, query_statistics_length: 'int | None' = None, max_amino_acid_distance: 'int | None' = None, mask_low_complexity: 'bool' = True, repeatmasker_dir: 'str | Path | None' = None, timeout_seconds: 'float' = 86400.0, verbose: 'bool' = True) -> 'dict[str, Any]'`。
- **输入：** 必填：`fasta`, `output_directory`；其余参数均为仅关键字参数且可选（`None` 或 `False` 表示使用工具自己的默认值）。
- **输出：** 含 `status, backend, operation, version, command, output_directory, outputs, log_tail, provenance` 的字典，另有该函数特有的键（见英文版 Output）。

## 14) `dnakit.annotation.repeatmasker_p11_cli_functions.mask_repeats()`

- **作用：** 筛查 DNA 中的散在重复和低复杂度序列，输出注释和掩蔽后的 FASTA（RepeatMasker）。
- **处理方法：** 运行固定的参数列表（不经过 shell）；参数含义、取值范围和工具默认值见英文版同名小节的 Method 说明。
- **API：** `dnakit.annotation.repeatmasker_p11_cli_functions.mask_repeats(fasta: 'str | Path', output_directory: 'str | Path', *, species: 'str | None' = None, library: 'str | Path | None' = None, search_engine: 'str | None' = None, speed: 'str | None' = None, parallel: 'int | None' = None, cutoff: 'int | None' = None, mask_low_complexity: 'bool' = True, mask_interspersed: 'bool' = True, screen_small_rnas: 'bool' = True, alu_only: 'bool' = False, max_divergence: 'float | None' = None, gc_level: 'int | None' = None, calculate_gc: 'bool' = False, max_fragment_length: 'int | None' = None, skip_excision: 'bool' = False, insertion_sequences: 'str' = 'report', mask_style: 'str' = 'n', lowercase_ambiguous: 'bool' = False, alignments: 'bool' = False, alignments_in_repeat_orientation: 'bool' = False, polymorphic_simple_repeats: 'bool' = False, ace_output: 'bool' = False, gff_output: 'bool' = False, html_output: 'bool' = False, untouched_annotation: 'bool' = False, crossmatch_output: 'bool' = False, no_id_column: 'bool' = False, exclude_n_runs: 'bool' = False, uncurated_families: 'bool' = False, cache_dir: 'str | Path | None' = None, famdb_data_dir: 'str | Path | None' = None, repeatmasker_dir: 'str | Path | None' = None, timeout_seconds: 'float' = 86400.0, verbose: 'bool' = True) -> 'dict[str, Any]'`。
- **输入：** 必填：`fasta`, `output_directory`；其余参数均为仅关键字参数且可选（`None` 或 `False` 表示使用工具自己的默认值）。
- **输出：** 含 `status, backend, operation, version, command, output_directory, outputs, log_tail, provenance` 的字典，另有该函数特有的键（见英文版 Output）。

## 15) `dnakit.annotation.repeatmasker_p11_cli_functions.mask_sequences_from_annotation()`

- **作用：** 根据重复注释对 FASTA 做硬掩蔽或软掩蔽（maskFile.pl）。
- **处理方法：** 运行固定的参数列表（不经过 shell）；参数含义、取值范围和工具默认值见英文版同名小节的 Method 说明。
- **API：** `dnakit.annotation.repeatmasker_p11_cli_functions.mask_sequences_from_annotation(fasta: 'str | Path', repeat_annotation: 'str | Path', output_directory: 'str | Path', *, soft_mask: 'bool' = False, min_divergence: 'float | None' = None, max_divergence: 'float | None' = None, repeatmasker_dir: 'str | Path | None' = None, timeout_seconds: 'float' = 3600.0, verbose: 'bool' = True) -> 'dict[str, Any]'`。
- **输入：** 必填：`fasta`, `repeat_annotation`, `output_directory`；其余参数均为仅关键字参数且可选（`None` 或 `False` 表示使用工具自己的默认值）。
- **输出：** 含 `status, backend, operation, version, command, output_directory, outputs, log_tail, provenance` 的字典，另有该函数特有的键（见英文版 Output）。

## 16) `dnakit.annotation.repeatmasker_p11_cli_functions.mask_tandem_repeats()`

- **作用：** 用 TRF 掩蔽 DNA 中的简单串联重复（trfMask）。
- **处理方法：** 运行固定的参数列表（不经过 shell）；参数含义、取值范围和工具默认值见英文版同名小节的 Method 说明。
- **API：** `dnakit.annotation.repeatmasker_p11_cli_functions.mask_tandem_repeats(fasta: 'str | Path', output_directory: 'str | Path', *, soft_mask: 'bool' = False, repeatmasker_dir: 'str | Path | None' = None, timeout_seconds: 'float' = 86400.0, verbose: 'bool' = True) -> 'dict[str, Any]'`。
- **输入：** 必填：`fasta`, `output_directory`；其余参数均为仅关键字参数且可选（`None` 或 `False` 表示使用工具自己的默认值）。
- **输出：** 含 `status, backend, operation, version, command, output_directory, outputs, log_tail, provenance` 的字典，另有该函数特有的键（见英文版 Output）。

## 17) `dnakit.annotation.repeatmasker_p11_cli_functions.plot_duplicon_map()`

- **作用：** 把 DupMasker 的 duplicon 图绘制为 SVG（dupliconToSVG.pl）。
- **处理方法：** 运行固定的参数列表（不经过 shell）；参数含义、取值范围和工具默认值见英文版同名小节的 Method 说明。
- **API：** `dnakit.annotation.repeatmasker_p11_cli_functions.plot_duplicon_map(duplicons_file: 'str | Path', output_directory: 'str | Path', *, repeatmasker_dir: 'str | Path | None' = None, timeout_seconds: 'float' = 3600.0, verbose: 'bool' = True) -> 'dict[str, Any]'`。
- **输入：** 必填：`duplicons_file`, `output_directory`；其余参数均为仅关键字参数且可选（`None` 或 `False` 表示使用工具自己的默认值）。
- **输出：** 含 `status, backend, operation, version, command, output_directory, outputs, log_tail, provenance` 的字典，另有该函数特有的键（见英文版 Output）。

## 18) `dnakit.annotation.repeatmasker_p11_cli_functions.plot_repeat_landscape()`

- **作用：** 绘制重复景观图：各重复类别在各分歧度上占基因组的比例（createRepeatLandscape.pl）。
- **处理方法：** 运行固定的参数列表（不经过 shell）；参数含义、取值范围和工具默认值见英文版同名小节的 Method 说明。
- **API：** `dnakit.annotation.repeatmasker_p11_cli_functions.plot_repeat_landscape(divergence_summary: 'str | Path', output_directory: 'str | Path', *, genome_size: 'int | None' = None, two_bit_file: 'str | Path | None' = None, title: 'str | None' = None, javascript_only: 'bool' = False, output_name: 'str' = 'landscape.html', repeatmasker_dir: 'str | Path | None' = None, timeout_seconds: 'float' = 3600.0, verbose: 'bool' = True) -> 'dict[str, Any]'`。
- **输入：** 必填：`divergence_summary`, `output_directory`；其余参数均为仅关键字参数且可选（`None` 或 `False` 表示使用工具自己的默认值）。
- **输出：** 含 `status, backend, operation, version, command, output_directory, outputs, log_tail, provenance` 的字典，另有该函数特有的键（见英文版 Output）。

## 19) `dnakit.annotation.repeatmasker_p11_cli_functions.process_repeat_hits()`

- **作用：** 把 RepeatMasker 原始搜索结果（.cat）后处理为重复注释（.out）和汇总（.tbl）。
- **处理方法：** 运行固定的参数列表（不经过 shell）；参数含义、取值范围和工具默认值见英文版同名小节的 Method 说明。
- **API：** `dnakit.annotation.repeatmasker_p11_cli_functions.process_repeat_hits(cat_file: 'str | Path', output_directory: 'str | Path', *, species: 'str | None' = None, library: 'str | Path | None' = None, mask_low_complexity: 'bool' = True, mask_interspersed: 'bool' = True, lowercase_ambiguous: 'bool' = False, untouched_annotation: 'bool' = False, crossmatch_output: 'bool' = False, ace_output: 'bool' = False, gff_output: 'bool' = False, polymorphic_simple_repeats: 'bool' = False, no_id_column: 'bool' = False, exclude_n_runs: 'bool' = False, align_l1_orf2: 'bool' = False, alignments: 'bool' = False, mask_source: 'str | Path | None' = None, mask_style: 'str' = 'n', repeatmasker_dir: 'str | Path | None' = None, timeout_seconds: 'float' = 86400.0, verbose: 'bool' = True) -> 'dict[str, Any]'`。
- **输入：** 必填：`cat_file`, `output_directory`；其余参数均为仅关键字参数且可选（`None` 或 `False` 表示使用工具自己的默认值）。
- **输出：** 含 `status, backend, operation, version, command, output_directory, outputs, log_tail, provenance` 的字典，另有该函数特有的键（见英文版 Output）。

## 20) `dnakit.annotation.repeatmasker_p11_cli_functions.recalculate_repeat_divergence()`

- **作用：** 重新计算重复比对的分歧度并按家族汇总（calcDivergenceFromAlign.pl）。
- **处理方法：** 运行固定的参数列表（不经过 shell）；参数含义、取值范围和工具默认值见英文版同名小节的 Method 说明。
- **API：** `dnakit.annotation.repeatmasker_p11_cli_functions.recalculate_repeat_divergence(alignment_file: 'str | Path', output_directory: 'str | Path', *, summary_name: 'str' = 'repeats.divsum', write_alignments: 'bool' = False, adjust_cpg: 'bool' = True, repeatmasker_dir: 'str | Path | None' = None, timeout_seconds: 'float' = 3600.0, verbose: 'bool' = True) -> 'dict[str, Any]'`。
- **输入：** 必填：`alignment_file`, `output_directory`；其余参数均为仅关键字参数且可选（`None` 或 `False` 表示使用工具自己的默认值）。
- **输出：** 含 `status, backend, operation, version, command, output_directory, outputs, log_tail, provenance` 的字典，另有该函数特有的键（见英文版 Output）。

## 21) `dnakit.annotation.repeatmasker_p11_cli_functions.summarize_repeat_annotation()`

- **作用：** 按家族和类别汇总重复注释（buildSummary.pl）。
- **处理方法：** 运行固定的参数列表（不经过 shell）；参数含义、取值范围和工具默认值见英文版同名小节的 Method 说明。
- **API：** `dnakit.annotation.repeatmasker_p11_cli_functions.summarize_repeat_annotation(repeat_annotation: 'str | Path', output_directory: 'str | Path', *, output_name: 'str' = 'summary.txt', min_divergence: 'float | None' = None, max_divergence: 'float | None' = None, species: 'str | None' = None, genome: 'str | Path | None' = None, use_absolute_genome_size: 'bool' = False, library_dir: 'str | Path | None' = None, repeatmasker_dir: 'str | Path | None' = None, timeout_seconds: 'float' = 3600.0, verbose: 'bool' = True) -> 'dict[str, Any]'`。
- **输入：** 必填：`repeat_annotation`, `output_directory`；其余参数均为仅关键字参数且可选（`None` 或 `False` 表示使用工具自己的默认值）。
- **输出：** 含 `status, backend, operation, version, command, output_directory, outputs, log_tail, provenance` 的字典，另有该函数特有的键（见英文版 Output）。

