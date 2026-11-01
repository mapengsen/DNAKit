# 序列文件与基因组注释

## 作用

为 P39/P40 整合批次提供明确的 DNA 操作接口，返回可序列化结果。

## 处理方法

函数调用其对应的上游操作；显式指定本地文件，新输出路径不得已存在。UCSC/黑名单下载仅在主动调用时联网。验证环境缺少 NUPACK 与 BLAST+，相关科学计算流程尚未实测。Enformer 使用与架构和适配器配置一致的本地权重。

## API

`dnakit.annotation.pyrodigal_p40_functions`

| Function | Operation |
| --- | --- |
| `pyrodigal_candidate_nodes` | Extract (and optionally score) the candidate start/stop nodes of a sequence. |
| `pyrodigal_codon_probability` | Start- or stop-codon probability expected from the GC content of a sequence. |
| `pyrodigal_find_genes` | Predict genes with pyrodigal and write the chosen annotation files. |
| `pyrodigal_gc_frame_plot` | Maximum-GC frame plot of one sequence (pyrodigal Sequence.max_gc_frame_plot). |
| `pyrodigal_shine_dalgarno` | Score the Shine-Dalgarno motif upstream of a start codon (pyrodigal Sequence.shine_dalgarno). |
| `pyrodigal_train_gene_finder` | Train the Prodigal gene finder on a genome FASTA (pyrodigal GeneFinder.train, ``pyrodigal -t``). |

`dnakit.sequence.pyfaidx_p40_functions`

| Function | Operation |
| --- | --- |
| `faidx_overwrite_region` | Overwrite a FASTA region with a sequence of the same length (pyfaidx Faidx.to_file). |
| `faidx_read_index` | Read an existing FASTA ``.fai`` index (pyfaidx Faidx.read_fai) or BGZF ``.gzi`` index (Faidx.read_gzi). |
| `faidx_write_index` | Build a FASTA index and write it to a new file. |

`dnakit.sequence.pyfastx_p40_functions`

| Function | Operation |
| --- | --- |
| `pyfastx_count_sequences` | Count the FASTA sequences that are at least ``min_length`` bp long (pyfastx Fasta.count). |
| `pyfastx_extract_records` | Extract FASTA sequences or FASTQ reads by name into a new file (``pyfastx extract``). |
| `pyfastx_list_sequence_names` | List the sequence names of a FASTA file (pyfastx Fasta.keys). |
| `pyfastx_nl_statistic` | Compute the N50/L50-style statistic of a FASTA file (pyfastx Fasta.nl). |
| `pyfastx_split_fastq` | Split a FASTQ file into several files (``pyfastx split``). |

`dnakit.sequence.dnaio_p40_functions`

| Function | Operation |
| --- | --- |
| `dnaio_read_parallel_files` | Read three or more parallel FASTA/FASTQ files record by record (dnaio.open with several files). |
| `dnaio_reverse_complement_reads` | Reverse-complement every read of a FASTA/FASTQ file (dnaio SequenceRecord.reverse_complement). |
| `dnaio_split_into_chunks` | Split a FASTA/FASTQ file into chunks of complete records, without parsing them (dnaio.read_chunks). |
| `dnaio_write_reads` | Write parallel read lists with dnaio: paired-end files, parallel FASTA/FASTQ files or interleaved reads. |

`dnakit.databases.genomepy_p40_functions`

| Function | Operation |
| --- | --- |
| `genomepy_list_installed` | List genomes installed under an explicit local directory, with no provider lookup. |
| `genomepy_import_local_annotation` | Copy a user-supplied local GTF/BED annotation to a new path and parse it with genomepy. |
| `genomepy_annotation_attributes` | List GTF attribute names using Annotation.attributes (upstream feature sampling). |
| `genomepy_annotation_genes` | List unique GTF gene_name values or BED names with Annotation.genes. |
| `genomepy_annotation_lengths` | Gene/transcript exon lengths with Annotation.lengths, preserving upstream coordinate arithmetic. |
| `genomepy_annotation_dictionary` | Make a key-value dictionary from GTF columns/attributes with Annotation.gtf_dict. |
| `genomepy_sanitize_annotation` | Match GTF contigs against full FASTA headers, then drop unmatched rows using genomepy sanitize. |
| `genomepy_download_blacklist` | Download and decompress the Boyle Lab blacklist for an assembly supported by genomepy (network). |

## 输入

DNA 字符串、JSON 字典/列表、本地基因组文件或模型权重。区间坐标、热力学参数见各函数 docstring；NUPACK 默认材料为 DNA。全部接口支持关键字参数 `verbose`。

## 输出

含完成状态、后端信息和操作结果的字典；文件操作返回输出路径，绘图返回图片路径及对应数据。缺少可选依赖时抛出 `BackendUnavailableError`。各批 manifest 单独记录未验收流程。

**示例代码：**

```python
from tempfile import TemporaryDirectory
from pathlib import Path
from dnakit.sequence.pyfastx_p40_functions import pyfastx_count_sequences
with TemporaryDirectory() as folder:
    fasta = Path(folder) / 'example.fa'
    fasta.write_text('>a\nACGT\n>b\nN\n')
    print(pyfastx_count_sequences(fasta, min_length=2, verbose=False)['sequence_count'])
```

**示例输出：**

```text
1
```

`pyrodigal_gc_frame_plot` 保留上游的未指定帧值 `-1`，在 `unknown_frame_count` 单独计数，并通过 `segments_truncated` 标记省略的区段。`dnaio_split_into_chunks` 在输出扩展名含压缩格式时会重新压缩每个分块。
