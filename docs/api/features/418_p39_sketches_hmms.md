# DNA 草图与 HMM

## 作用

为 P39/P40 整合批次提供明确的 DNA 操作接口，返回可序列化结果。

## 处理方法

函数调用其对应的上游操作；显式指定本地文件，新输出路径不得已存在。UCSC/黑名单下载仅在主动调用时联网。验证环境缺少 NUPACK 与 BLAST+，相关科学计算流程尚未实测。Enformer 使用与架构和适配器配置一致的本地权重。

## API

`dnakit.alignment.pyhmmer_p39_functions`

| Function | Operation |
| --- | --- |
| `build_hmm_from_sequence` | Build a DNA profile HMM from one sequence and write it as a HMMER3 file (pyhmmer.plan7.Builder.build). |
| `build_sequence_ssi_index` | Index the records of a FASTA file by name in a new SSI file next to it (SSIWriter.add_file, add_key, add_alias). |
| `compute_alignment_weights` | Compute sequence weights of a DNA alignment with Easel (MSA.compute_weights, esl-weight). |
| `emit_hmm_alignment` | Sample a multiple alignment of sequences emitted from the core model of a DNA profile HMM (HMM.emit_alignment). |
| `emit_hmm_sequences` | Sample sequences from a DNA profile HMM and write them to a new FASTA file (HMM.emit_sequence, Profile.emit_sequence). |
| `fetch_sequence_by_ssi` | Fetch one record of a FASTA file by name or alias through its SSI index (SSIReader.find_name). |
| `filter_alignment_by_identity` | Remove alignment sequences whose pairwise identity reaches ``max_identity`` and write the rest to a new file (DigitalMSA.identity_filter, esl-weight -f). |
| `guess_file_alphabet` | Detect whether a sequence file or an alignment file holds DNA, RNA or amino acids (SequenceFile.guess_alphabet, MSAFile.guess_alphabet). |
| `pyhmmer_hmmscan` | Scan each DNA sequence of a FASTA file against a database of DNA profile HMMs with the HMMER 3 hmmscan pipeline (pyhmmer.hmmscan). |
| `pyhmmer_hmmsearch` | Search a DNA sequence file with one or more DNA profile HMMs using the HMMER 3 hmmsearch pipeline (pyhmmer.hmmsearch). |
| `reverse_complement_alignment` | Reverse-complement a DNA multiple alignment, including per-column annotation, and write it to a new file (DigitalMSA.reverse_complement). |

`dnakit.microbiome.sourmash_p39_functions`

| Function | Operation |
| --- | --- |
| `build_signature_index` | Build an on-disk searchable signature index of a chosen type (sourmash index -F). |
| `ingest_mash_json_signatures` | Convert JSON sketch dumps written by ``mash info -d`` into one sourmash signature file (sourmash sig ingest). |
| `murmur_hash_string` | MurmurHash64 of a string exactly as given, with sourmash's seed 42 by default (sourmash.minhash.hash_murmur). |
| `sequence_kmer_hashes` | Hash every k-mer of a DNA sequence without making a sketch (MinHash.seq_to_hashes, kmers_and_hashes). |
| `sketch_intersection_union_size` | Number of shared hashes and the size of the union of two sketches (MinHash.count_common and intersection_and_union_size). |
| `sketch_size_accuracy` | Tell whether a scaled sketch is large enough to estimate its k-mer set size (MinHash.size_is_accurate). |
| `sketch_weighted_containment` | Abundance-weighted fraction of the query sketch found in the target sketch (MinHash.contained_by_weighted). |
| `tax_annotate_gather_results` | Add a lineage column to gather CSV files (sourmash tax annotate); writes ``<name>.with-lineages.csv``. |
| `tax_genome_classification` | Assign one lineage to each query genome from gather results (sourmash tax genome). |
| `tax_grep_lineages` | Select the taxonomy lineages whose rank names match a string (sourmash tax grep). |
| `tax_prepare_taxonomy_database` | Combine taxonomy spreadsheets into one new taxonomy database, SQLite (``sql``) or CSV (sourmash tax prepare). |
| `tax_summarize_lineages` | Summarize taxonomy databases or lineage spreadsheets: number of distinct lineages per rank (sourmash tax summarize). |
| `watch_reads_against_index` | Stream FASTA/FASTQ reads into a growing sketch and stop at the first index match (sourmash watch). |

## 输入

DNA 字符串、JSON 字典/列表、本地基因组文件或模型权重。区间坐标、热力学参数见各函数 docstring；NUPACK 默认材料为 DNA。全部接口支持关键字参数 `verbose`。

## 输出

含完成状态、后端信息和操作结果的字典；文件操作返回输出路径，绘图返回图片路径及对应数据。缺少可选依赖时抛出 `BackendUnavailableError`。各批 manifest 单独记录未验收流程。

**示例代码：**

```python
from dnakit.microbiome.sourmash_p39_functions import murmur_hash_string
print(murmur_hash_string('ACGT', verbose=False)['hash'])
```

**示例输出：**

```text
2597925387403686983
```
