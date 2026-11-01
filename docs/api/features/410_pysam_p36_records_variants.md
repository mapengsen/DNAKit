# pysam 比对与变异记录

本页记录 P36 相关的真实 DNAKit 接口。函数按所需可选依赖延迟导入；坐标、默认值与参数选择以各接口说明为准。

验证范围说明：VCF positions are 1-based; BAM intervals use 0-based half-open coordinates. Legacy samtools pileup -c text conversion is excluded.

## `dnakit.alignment.pysam_p36_alignment_functions.fetch_mate_alignments`

**作用：** Find the mate alignment of every paired, mapped read in a region (``AlignmentFile.mate``).

**处理方法：** 由模块中的显式接口校验输入并调用对应后端。

**API：**

```text
def fetch_mate_alignments(alignment_file: str | Path, contig: str, start: int | None=None, end: int | None=None, *, max_records: int=1000, verbose: bool=True) -> dict[str, Any]:
```

**输入：** 签名列出必选参数、可选参数与原生默认值；参数含义和输出如下。

```text
Find the mate alignment of every paired, mapped read in a region (``AlignmentFile.mate``).

The file must be indexed and coordinate sorted. Reads that are unpaired or whose mate is unmapped are
listed in ``without_mate``; every other read gets the mate's ``pysam`` record dictionary.
```

**输出：** 返回可 JSON 序列化的字典；文件输出返回路径。可选依赖缺失时抛出带安装提示的 BackendUnavailableError。

## `dnakit.alignment.pysam_p36_alignment_functions.count_base_coverage`

**作用：** Count A, C, G and T bases per position of a region (``AlignmentFile.count_coverage``).

**处理方法：** 由模块中的显式接口校验输入并调用对应后端。

**API：**

```text
def count_base_coverage(alignment_file: str | Path, contig: str, start: int, end: int, *, quality_threshold: int=15, read_callback: str='all', verbose: bool=True) -> dict[str, Any]:
```

**输入：** 签名列出必选参数、可选参数与原生默认值；参数含义和输出如下。

```text
Count A, C, G and T bases per position of a region (``AlignmentFile.count_coverage``).

``quality_threshold`` skips bases below this Phred quality (pysam default 15); ``read_callback`` is
``all`` (skip unmapped, secondary, QC-fail and duplicate reads) or ``nofilter``. The four lists have one
count per position from ``start`` to ``end``.
```

**输出：** 返回可 JSON 序列化的字典；文件输出返回路径。可选依赖缺失时抛出带安装提示的 BackendUnavailableError。

## `dnakit.alignment.pysam_p36_alignment_functions.count_spliced_junctions`

**作用：** Count reads that support each intron (CIGAR ``N`` gap) in a region (``AlignmentFile.find_introns``).

**处理方法：** 由模块中的显式接口校验输入并调用对应后端。

**API：**

```text
def count_spliced_junctions(alignment_file: str | Path, contig: str | None=None, start: int | None=None, end: int | None=None, *, method: str='fast', verbose: bool=True) -> dict[str, Any]:
```

**输入：** 签名列出必选参数、可选参数与原生默认值；参数含义和输出如下。

```text
Count reads that support each intron (CIGAR ``N`` gap) in a region (``AlignmentFile.find_introns``).

``method`` is ``fast`` (``find_introns``) or ``slow`` (``find_introns_slow``, which collects the
reads' aligned blocks); both give the same junctions. Each junction is the 0-based half-open
interval of the intron with its read count.
```

**输出：** 返回可 JSON 序列化的字典；文件输出返回路径。可选依赖缺失时抛出带安装提示的 BackendUnavailableError。

## `dnakit.alignment.pysam_p36_alignment_functions.pileup_columns`

**作用：** Iterate the pileup columns of a region and report the reads aligned at each (``AlignmentFile.pileup``).

**处理方法：** 由模块中的显式接口校验输入并调用对应后端。

**API：**

```text
def pileup_columns(alignment_file: str | Path, contig: str, start: int | None=None, end: int | None=None, *, min_base_quality: int=13, min_mapping_quality: int=0, stepper: str='all', truncate: bool=False, max_depth: int=8000, ignore_overlaps: bool=True, ignore_orphans: bool=True, include_read_names: bool=False, max_columns: int=100000, verbose: bool=True) -> dict[str, Any]:
```

**输入：** 签名列出必选参数、可选参数与原生默认值；参数含义和输出如下。

```text
Iterate the pileup columns of a region and report the reads aligned at each (``AlignmentFile.pileup``).

For every column: ``reference_pos`` (0-based), ``nsegments`` (reads covering it, before the base-quality
filter), ``num_aligned`` (``PileupColumn.get_num_aligned`` - reads whose base passes
``min_base_quality``), and optionally the read names. ``stepper`` is ``all``, ``nofilter``, ``samtools``
or ``bam``; ``truncate=False`` is pysam's default and may list columns just outside the region.
```

**输出：** 返回可 JSON 序列化的字典；文件输出返回路径。可选依赖缺失时抛出带安装提示的 BackendUnavailableError。

## `dnakit.alignment.pysam_p36_alignment_functions.find_alignments_by_name`

**作用：** Look up all alignments of the given read names with pysam's in-memory name index (``IndexedReads``).

**处理方法：** 由模块中的显式接口校验输入并调用对应后端。

**API：**

```text
def find_alignments_by_name(alignment_file: str | Path, query_names: Sequence[str], *, verbose: bool=True) -> dict[str, Any]:
```

**输入：** 签名列出必选参数、可选参数与原生默认值；参数含义和输出如下。

```text
Look up all alignments of the given read names with pysam's in-memory name index (``IndexedReads``).

The whole file is read once to build the index, so it needs no BAM index and works on unsorted files.
Names without an alignment map to an empty list.
```

**输出：** 返回可 JSON 序列化的字典；文件输出返回路径。可选依赖缺失时抛出带安装提示的 BackendUnavailableError。

## `dnakit.alignment.pysam_p36_alignment_functions.sam_header_to_dict`

**作用：** Return the header of a SAM, BAM or CRAM file as nested dictionaries (``AlignmentHeader.to_dict``).

**处理方法：** 由模块中的显式接口校验输入并调用对应后端。

**API：**

```text
def sam_header_to_dict(alignment_file: str | Path, *, verbose: bool=True) -> dict[str, Any]:
```

**输入：** 签名列出必选参数、可选参数与原生默认值；参数含义和输出如下。

```text
Return the header of a SAM, BAM or CRAM file as nested dictionaries (``AlignmentHeader.to_dict``).
```

**输出：** 返回可 JSON 序列化的字典；文件输出返回路径。可选依赖缺失时抛出带安装提示的 BackendUnavailableError。

## `dnakit.alignment.pysam_p36_alignment_functions.write_sam_header`

**作用：** Build a SAM header and write it as text (``AlignmentHeader.from_dict`` / ``from_references``).

**处理方法：** 由模块中的显式接口校验输入并调用对应后端。

**API：**

```text
def write_sam_header(output_path: str | Path, *, header_dict: Mapping[str, Any] | None=None, reference_names: Sequence[str] | None=None, reference_lengths: Sequence[int] | None=None, verbose: bool=True) -> dict[str, Any]:
```

**输入：** 签名列出必选参数、可选参数与原生默认值；参数含义和输出如下。

```text
Build a SAM header and write it as text (``AlignmentHeader.from_dict`` / ``from_references``).

Give either ``header_dict`` (the nested form of ``sam_header_to_dict``: ``HD``, ``SQ`` list, ``RG`` list,
``PG`` list, ``CO`` list) or ``reference_names`` with ``reference_lengths`` (one ``@SQ`` line each).
The result holds the header text and the round-tripped dictionary.
```

**输出：** 返回可 JSON 序列化的字典；文件输出返回路径。可选依赖缺失时抛出带安装提示的 BackendUnavailableError。

## `dnakit.alignment.pysam_p36_alignment_functions.alignment_records_to_dicts`

**作用：** Convert alignment records to SAM-style dictionaries (``AlignedSegment.to_dict``).

**处理方法：** 由模块中的显式接口校验输入并调用对应后端。

**API：**

```text
def alignment_records_to_dicts(alignment_file: str | Path, contig: str | None=None, start: int | None=None, end: int | None=None, *, max_records: int=1000, verbose: bool=True) -> dict[str, Any]:
```

**输入：** 签名列出必选参数、可选参数与原生默认值；参数含义和输出如下。

```text
Convert alignment records to SAM-style dictionaries (``AlignedSegment.to_dict``).

Each record holds ``name, flag, ref_name, ref_pos`` (1-based), ``map_quality, cigar, next_ref_name,
next_ref_pos, length, seq, qual`` and ``tags`` - all as SAM text, so the output can be written back with
``write_alignments_from_dicts``.
```

**输出：** 返回可 JSON 序列化的字典；文件输出返回路径。可选依赖缺失时抛出带安装提示的 BackendUnavailableError。

## `dnakit.alignment.pysam_p36_alignment_functions.write_alignments_from_dicts`

**作用：** Build alignment records from SAM-style dictionaries and write them (``AlignedSegment.from_dict``).

**处理方法：** 由模块中的显式接口校验输入并调用对应后端。

**API：**

```text
def write_alignments_from_dicts(records: Sequence[Mapping[str, Any]], header: Mapping[str, Any], output_path: str | Path, *, verbose: bool=True) -> dict[str, Any]:
```

**输入：** 签名列出必选参数、可选参数与原生默认值；参数含义和输出如下。

```text
Build alignment records from SAM-style dictionaries and write them (``AlignedSegment.from_dict``).

``header`` is the nested header dictionary; ``records`` use the keys of ``alignment_records_to_dicts``.
The output format follows the suffix (``.sam`` or ``.bam``).
```

**输出：** 返回可 JSON 序列化的字典；文件输出返回路径。可选依赖缺失时抛出带安装提示的 BackendUnavailableError。

## `dnakit.alignment.pysam_p36_alignment_functions.read_aligned_pairs`

**作用：** List the (read position, reference position) pairs of each read (``AlignedSegment.get_aligned_pairs``).

**处理方法：** 由模块中的显式接口校验输入并调用对应后端。

**API：**

```text
def read_aligned_pairs(alignment_file: str | Path, contig: str | None=None, start: int | None=None, end: int | None=None, *, matches_only: bool=False, with_seq: bool=False, max_records: int=1000, verbose: bool=True) -> dict[str, Any]:
```

**输入：** 签名列出必选参数、可选参数与原生默认值；参数含义和输出如下。

```text
List the (read position, reference position) pairs of each read (``AlignedSegment.get_aligned_pairs``).

Insertions and soft clips have ``None`` as the reference position, deletions and skipped regions ``None``
as the read position, unless ``matches_only=True``. ``with_seq=True`` adds the reference base (needs an MD
tag) as third element.
```

**输出：** 返回可 JSON 序列化的字典；文件输出返回路径。可选依赖缺失时抛出带安装提示的 BackendUnavailableError。

## `dnakit.alignment.pysam_p36_alignment_functions.read_aligned_blocks`

**作用：** List the gapless aligned blocks of each read in reference coordinates (``AlignedSegment.get_blocks``).

**处理方法：** 由模块中的显式接口校验输入并调用对应后端。

**API：**

```text
def read_aligned_blocks(alignment_file: str | Path, contig: str | None=None, start: int | None=None, end: int | None=None, *, max_records: int=1000, verbose: bool=True) -> dict[str, Any]:
```

**输入：** 签名列出必选参数、可选参数与原生默认值；参数含义和输出如下。

```text
List the gapless aligned blocks of each read in reference coordinates (``AlignedSegment.get_blocks``).

A block is a 0-based half-open ``[start, end)`` interval; deletions and introns split blocks.
```

**输出：** 返回可 JSON 序列化的字典；文件输出返回路径。可选依赖缺失时抛出带安装提示的 BackendUnavailableError。

## `dnakit.alignment.pysam_p36_alignment_functions.read_cigar_stats`

**作用：** Summarise the CIGAR of each read: bases and blocks per operation, plus the NM edit distance

**处理方法：** 由模块中的显式接口校验输入并调用对应后端。

**API：**

```text
def read_cigar_stats(alignment_file: str | Path, contig: str | None=None, start: int | None=None, end: int | None=None, *, max_records: int=1000, verbose: bool=True) -> dict[str, Any]:
```

**输入：** 签名列出必选参数、可选参数与原生默认值；参数含义和输出如下。

```text
Summarise the CIGAR of each read: bases and blocks per operation, plus the NM edit distance
(``AlignedSegment.get_cigar_stats``).

``bases`` and ``blocks`` map the operation names M, I, D, N, S, H, P, =, X and NM to counts.
```

**输出：** 返回可 JSON 序列化的字典；文件输出返回路径。可选依赖缺失时抛出带安装提示的 BackendUnavailableError。

## `dnakit.alignment.pysam_p36_alignment_functions.read_reference_overlap`

**作用：** Count how many aligned bases of each read overlap ``[start, end)`` on the reference

**处理方法：** 由模块中的显式接口校验输入并调用对应后端。

**API：**

```text
def read_reference_overlap(alignment_file: str | Path, contig: str, start: int, end: int, *, max_records: int=1000, verbose: bool=True) -> dict[str, Any]:
```

**输入：** 签名列出必选参数、可选参数与原生默认值；参数含义和输出如下。

```text
Count how many aligned bases of each read overlap ``[start, end)`` on the reference
(``AlignedSegment.get_overlap``).
```

**输出：** 返回可 JSON 序列化的字典；文件输出返回路径。可选依赖缺失时抛出带安装提示的 BackendUnavailableError。

## `dnakit.alignment.pysam_p36_alignment_functions.read_reference_positions`

**作用：** List the 0-based reference positions covered by aligned read bases (``AlignedSegment.get_reference_positions``).

**处理方法：** 由模块中的显式接口校验输入并调用对应后端。

**API：**

```text
def read_reference_positions(alignment_file: str | Path, contig: str | None=None, start: int | None=None, end: int | None=None, *, full_length: bool=False, max_records: int=1000, verbose: bool=True) -> dict[str, Any]:
```

**输入：** 签名列出必选参数、可选参数与原生默认值；参数含义和输出如下。

```text
List the 0-based reference positions covered by aligned read bases (``AlignedSegment.get_reference_positions``).

``full_length=True`` returns one entry per read base with ``None`` for insertions and clips.
```

**输出：** 返回可 JSON 序列化的字典；文件输出返回路径。可选依赖缺失时抛出带安装提示的 BackendUnavailableError。

## `dnakit.alignment.pysam_p36_alignment_functions.reconstruct_read_reference_sequence`

**作用：** Rebuild the reference sequence under each read from its MD tag (``AlignedSegment.get_reference_sequence``).

**处理方法：** 由模块中的显式接口校验输入并调用对应后端。

**API：**

```text
def reconstruct_read_reference_sequence(alignment_file: str | Path, contig: str | None=None, start: int | None=None, end: int | None=None, *, max_records: int=1000, verbose: bool=True) -> dict[str, Any]:
```

**输入：** 签名列出必选参数、可选参数与原生默认值；参数含义和输出如下。

```text
Rebuild the reference sequence under each read from its MD tag (``AlignedSegment.get_reference_sequence``).

Reference bases at mismatches are lower case, matches upper case, as pysam returns them. Reads without an MD tag are
listed in ``skipped_without_md`` instead of failing.
```

**输出：** 返回可 JSON 序列化的字典；文件输出返回路径。可选依赖缺失时抛出带安装提示的 BackendUnavailableError。

## `dnakit.alignment.pysam_p36_alignment_functions.infer_read_lengths`

**作用：** Infer each read's length from its CIGAR (``infer_query_length`` and ``infer_read_length``).

**处理方法：** 由模块中的显式接口校验输入并调用对应后端。

**API：**

```text
def infer_read_lengths(alignment_file: str | Path, contig: str | None=None, start: int | None=None, end: int | None=None, *, max_records: int=1000, verbose: bool=True) -> dict[str, Any]:
```

**输入：** 签名列出必选参数、可选参数与原生默认值；参数含义和输出如下。

```text
Infer each read's length from its CIGAR (``infer_query_length`` and ``infer_read_length``).

``query_length`` counts M, I, S, =, X (without hard clips); ``read_length`` also counts hard-clipped
bases (the original read length). ``stored_length`` is the length of the stored sequence (None if
the sequence is ``*``).
```

**输出：** 返回可 JSON 序列化的字典；文件输出返回路径。可选依赖缺失时抛出带安装提示的 BackendUnavailableError。

## `dnakit.alignment.pysam_p36_alignment_functions.read_base_modifications`

**作用：** Decode the base modifications of each read from its MM/ML tags (``modified_bases``).

**处理方法：** 由模块中的显式接口校验输入并调用对应后端。

**API：**

```text
def read_base_modifications(alignment_file: str | Path, contig: str | None=None, start: int | None=None, end: int | None=None, *, forward: bool=False, max_records: int=1000, verbose: bool=True) -> dict[str, Any]:
```

**输入：** 签名列出必选参数、可选参数与原生默认值；参数含义和输出如下。

```text
Decode the base modifications of each read from its MM/ML tags (``modified_bases``).

For every modification (original base, strand, code such as ``m``) the result lists
``[read position, probability 0-255]`` pairs. ``forward=True`` uses ``modified_bases_forward`` (positions
in the original sequencing direction). Reads without MM/ML tags have an empty mapping.
```

**输出：** 返回可 JSON 序列化的字典；文件输出返回路径。可选依赖缺失时抛出带安装提示的 BackendUnavailableError。

## `dnakit.variants.pysam_p36_variant_functions.write_vcf_records`

**作用：** Build a VCF header and write VCF, bgzipped VCF or BCF records (``VariantFile.write``).

**处理方法：** 由模块中的显式接口校验输入并调用对应后端。

**API：**

```text
def write_vcf_records(output_path: str | Path, records: Sequence[Mapping[str, Any]], *, contigs: Sequence[Mapping[str, Any]]=(), info_fields: Sequence[Mapping[str, Any]]=(), format_fields: Sequence[Mapping[str, Any]]=(), filters: Sequence[Mapping[str, Any]]=(), meta: Sequence[Mapping[str, Any]]=(), header_lines: Sequence[str]=(), samples: Sequence[str]=(), verbose: bool=True) -> dict[str, Any]:
```

**输入：** 签名列出必选参数、可选参数与原生默认值；参数含义和输出如下。

```text
Build a VCF header and write VCF, bgzipped VCF or BCF records (``VariantFile.write``).

The header is made with ``VariantHeader.add_meta`` / ``add_line`` / ``add_sample``: ``contigs`` are
``{"id", "length"}``; ``info_fields`` and ``format_fields`` are ``{"id", "number", "type",
"description"}``; ``filters`` are ``{"id", "description"}``; ``meta`` are ``{"key", "value"}`` lines such
as ``fileDate``; ``header_lines`` are raw ``##...`` lines; ``samples`` are the sample names.
Each record is ``{"contig", "pos", "ref", "alts", "id", "qual", "filter", "info", "samples"}`` with a
1-based ``pos``; ``samples`` maps a sample name to its FORMAT values (``GT`` as a list such as
``[0, 1]``; ``"phased": true`` marks a phased genotype). The output suffix selects the format.
```

**输出：** 返回可 JSON 序列化的字典；文件输出返回路径。可选依赖缺失时抛出带安装提示的 BackendUnavailableError。

## `dnakit.variants.pysam_p36_variant_functions.merge_vcf_headers`

**作用：** Merge the header of a second VCF into a copy of the first header and write it (``VariantHeader.merge``).

**处理方法：** 由模块中的显式接口校验输入并调用对应后端。

**API：**

```text
def merge_vcf_headers(first_vcf: str | Path, second_vcf: str | Path, output_path: str | Path, *, verbose: bool=True) -> dict[str, Any]:
```

**输入：** 签名列出必选参数、可选参数与原生默认值；参数含义和输出如下。

```text
Merge the header of a second VCF into a copy of the first header and write it (``VariantHeader.merge``).

The output is a header-only VCF. Lines present in the second header but not in the first are added;
pysam does not merge sample columns, so the samples of the first header are kept.
```

**输出：** 返回可 JSON 序列化的字典；文件输出返回路径。可选依赖缺失时抛出带安装提示的 BackendUnavailableError。

## `dnakit.variants.pysam_p36_variant_functions.translate_vcf_records`

**作用：** Rewrite the records of one VCF under the header of another VCF (``VariantRecord.translate``).

**处理方法：** 由模块中的显式接口校验输入并调用对应后端。

**API：**

```text
def translate_vcf_records(source_vcf: str | Path, header_vcf: str | Path, output_path: str | Path, *, max_records: int=100000, verbose: bool=True) -> dict[str, Any]:
```

**输入：** 签名列出必选参数、可选参数与原生默认值；参数含义和输出如下。

```text
Rewrite the records of one VCF under the header of another VCF (``VariantRecord.translate``).

Contigs, INFO, FORMAT and FILTER definitions are looked up in the target header; a record that uses a
definition missing there fails with an error that names it. Sample columns of the target header that
the source lacks stay missing.
```

**输出：** 返回可 JSON 序列化的字典；文件输出返回路径。可选依赖缺失时抛出带安装提示的 BackendUnavailableError。

## `dnakit.variants.pysam_p36_variant_functions.compare_variant_calls`

**作用：** Decide whether two variant calls are the same event (``pysam.VCF.compare_calls``).

**处理方法：** 由模块中的显式接口校验输入并调用对应后端。

**API：**

```text
def compare_variant_calls(pos1: int, ref1: str, alt1: str, pos2: int, ref2: str, alt2: str, *, verbose: bool=True) -> dict[str, Any]:
```

**输入：** 签名列出必选参数、可选参数与原生默认值；参数含义和输出如下。

```text
Decide whether two variant calls are the same event (``pysam.VCF.compare_calls``).

The positions must be equal; then the shared rightmost bases of REF and ALT are trimmed from each call
and only the remaining ALT strings are compared (pysam does not compare the trimmed REF strings).
``ATT>AT`` and ``AT>A`` at the same position are the same deletion of one T; ``A>C`` and ``A>G`` are not.
```

**输出：** 返回可 JSON 序列化的字典；文件输出返回路径。可选依赖缺失时抛出带安装提示的 BackendUnavailableError。

## `dnakit.variants.pysam_p36_variant_functions.convert_iupac_genotype_codes`

**作用：** Convert diploid genotypes to one-letter IUPAC codes and back (``Pileup.encodeGenotype`` / ``decodeGenotype``).

**处理方法：** 由模块中的显式接口校验输入并调用对应后端。

**API：**

```text
def convert_iupac_genotype_codes(values: Sequence[str], *, direction: str='encode', verbose: bool=True) -> dict[str, Any]:
```

**输入：** 签名列出必选参数、可选参数与原生默认值；参数含义和输出如下。

```text
Convert diploid genotypes to one-letter IUPAC codes and back (``Pileup.encodeGenotype`` / ``decodeGenotype``).

``direction="encode"`` maps two-base genotypes such as ``AG`` to ``r`` (lower case when the first base
sorts before the second, upper case otherwise; homozygous ``GG`` gives ``G``). ``direction="decode"``
maps codes such as ``m`` or ``M`` back to the two-base genotype. Unknown values are reported in
``errors`` and give ``None``.
```

**输出：** 返回可 JSON 序列化的字典；文件输出返回路径。可选依赖缺失时抛出带安装提示的 BackendUnavailableError。

## `dnakit.variants.pysam_p36_variant_functions.convert_gff_attributes`

**作用：** Convert the attribute column of a GTF or GFF3 line between its text form and a dictionary.

**处理方法：** 由模块中的显式接口校验输入并调用对应后端。

**API：**

```text
def convert_gff_attributes(*, attribute_string: str | None=None, attributes: Mapping[str, Any] | None=None, flavor: str='gtf', verbose: bool=True) -> dict[str, Any]:
```

**输入：** 签名列出必选参数、可选参数与原生默认值；参数含义和输出如下。

```text
Convert the attribute column of a GTF or GFF3 line between its text form and a dictionary.

Give ``attribute_string`` (for example ``gene_id "A"; gene_name "B";`` for ``flavor="gtf"`` or
``ID=a;Name=b`` for ``flavor="gff3"``) to get the parsed dictionary (``GTFProxy.attribute_string2dict`` /
``GFF3Proxy.attribute_string2dict``), or ``attributes`` to get the formatted string
(``dict2attribute_string``). GFF3 values stay percent-encoded as in the file; the string form does not
add escapes.
```

**输出：** 返回可 JSON 序列化的字典；文件输出返回路径。可选依赖缺失时抛出带安装提示的 BackendUnavailableError。

## 可执行示例

**示例代码：**

```python
from dnakit.variants.pysam_p36_variant_functions import convert_iupac_genotype_codes
print(convert_iupac_genotype_codes(['AG', 'GG'], verbose=False)['results'])
```

**示例输出：**

```text
{'AG': 'r', 'GG': 'G'}
```
