# PyEnsembl 本地注释查询

本页记录 P36 相关的真实 DNAKit 接口。函数按所需可选依赖延迟导入；坐标、默认值与参数选择以各接口说明为准。

## `dnakit.annotation.pyensembl_p36_functions.load_hgnc_gene_aliases`

**作用：** Load the HGNC gene-name alias table (symbols, aliases, previous symbols) mapped to Ensembl gene IDs.

**处理方法：** 由模块中的显式接口校验输入并调用对应后端。

**API：**

```text
def load_hgnc_gene_aliases(hgnc_tsv: str | Path | None=None, *, download: bool=False, source_url: str | None=None, cache_directory: str | Path | None=None, names: Sequence[str] | None=None, max_names: int=100000, verbose: bool=True) -> dict[str, Any]:
```

**输入：** 签名列出必选参数、可选参数与原生默认值；参数含义和输出如下。

```text
Load the HGNC gene-name alias table (symbols, aliases, previous symbols) mapped to Ensembl gene IDs.

With ``hgnc_tsv`` the file is read offline (``GeneNameAliases.from_hgnc``; plain or gzip complete-set
TSV). With ``download=True`` the current HGNC complete set (or ``source_url``) is fetched and cached
(``GeneNameAliases.download_hgnc``; network). ``names`` restricts the result to those symbols; otherwise
at most ``max_names`` entries are returned and ``total_names`` gives the full size.
```

**输出：** 返回可 JSON 序列化的字典；文件输出返回路径。可选依赖缺失时抛出带安装提示的 BackendUnavailableError。

## `dnakit.annotation.pyensembl_p36_functions.fetch_genomic_sequence`

**作用：** Read genomic DNA for a 1-based inclusive interval from a local genome FASTA (``Genome.sequence``).

**处理方法：** 由模块中的显式接口校验输入并调用对应后端。

**API：**

```text
def fetch_genomic_sequence(genome_fasta: str | Path, contig: str, start: int, end: int, *, strand: str='+', mask: str='upper', verbose: bool=True) -> dict[str, Any]:
```

**输入：** 签名列出必选参数、可选参数与原生默认值；参数含义和输出如下。

```text
Read genomic DNA for a 1-based inclusive interval from a local genome FASTA (``Genome.sequence``).

``strand="-"`` returns the reverse complement. ``mask="upper"`` upper-cases the result, ``mask="raw"``
keeps soft masking (pyensembl's default is ``upper``). The contig may be named as in the FASTA.
```

**输出：** 返回可 JSON 序列化的字典；文件输出返回路径。可选依赖缺失时抛出带安装提示的 BackendUnavailableError。

## `dnakit.annotation.pyensembl_p36_functions.fetch_transcript_sequence`

**作用：** Look up the cDNA sequences of transcript IDs in a local transcript FASTA (``Genome.transcript_sequence``).

**处理方法：** 由模块中的显式接口校验输入并调用对应后端。

**API：**

```text
def fetch_transcript_sequence(transcript_fasta: str | Path, transcript_ids: Sequence[str], *, verbose: bool=True) -> dict[str, Any]:
```

**输入：** 签名列出必选参数、可选参数与原生默认值；参数含义和输出如下。

```text
Look up the cDNA sequences of transcript IDs in a local transcript FASTA (``Genome.transcript_sequence``).

IDs absent from the FASTA map to ``null`` in ``sequences`` and are listed in ``missing``. A versioned ID
such as ``ENST00000269305.8`` is accepted when the FASTA records that version.
```

**输出：** 返回可 JSON 序列化的字典；文件输出返回路径。可选依赖缺失时抛出带安装提示的 BackendUnavailableError。

## `dnakit.annotation.pyensembl_p36_functions.query_annotation_at_locus`

**作用：** List the transcripts, exons or protein IDs that overlap a locus (``Genome.transcripts_at_locus``,

**处理方法：** 由模块中的显式接口校验输入并调用对应后端。

**API：**

```text
def query_annotation_at_locus(gtf: str | Path, contig: str, start: int, end: int | None=None, *, feature: str='transcript', strand: str | None=None, output: str='records', cache_dir: str | Path | None=None, max_results: int=1000, verbose: bool=True) -> dict[str, Any]:
```

**输入：** 签名列出必选参数、可选参数与原生默认值；参数含义和输出如下。

```text
List the transcripts, exons or protein IDs that overlap a locus (``Genome.transcripts_at_locus``,
``exons_at_locus``, ``protein_ids_at_locus`` and their ``*_ids_at_locus`` / ``*_names_at_locus`` forms).

``feature`` is ``transcript``, ``exon`` or ``protein``; ``output`` is ``records`` (objects as
dictionaries), ``ids``, or ``names`` (transcripts only). ``end`` defaults to ``start``. For genes use the
existing ``query_genes_at_locus``.
```

**输出：** 返回可 JSON 序列化的字典；文件输出返回路径。可选依赖缺失时抛出带安装提示的 BackendUnavailableError。

## `dnakit.annotation.pyensembl_p36_functions.merge_gene_intervals`

**作用：** Union of all gene loci on a contig as sorted non-overlapping intervals (``Genome.merged_gene_intervals``).

**处理方法：** 由模块中的显式接口校验输入并调用对应后端。

**API：**

```text
def merge_gene_intervals(gtf: str | Path, contig: str, *, strand: str | None=None, cache_dir: str | Path | None=None, verbose: bool=True) -> dict[str, Any]:
```

**输入：** 签名列出必选参数、可选参数与原生默认值；参数含义和输出如下。

```text
Union of all gene loci on a contig as sorted non-overlapping intervals (``Genome.merged_gene_intervals``).

Adjacent intervals (``end + 1 == next start``) are merged; ``strand`` restricts to one strand.
```

**输出：** 返回可 JSON 序列化的字典；文件输出返回路径。可选依赖缺失时抛出带安装提示的 BackendUnavailableError。

## `dnakit.annotation.pyensembl_p36_functions.find_nearest_transcript`

**作用：** Find the transcript whose locus is nearest to a position or interval, even if none overlaps

**处理方法：** 由模块中的显式接口校验输入并调用对应后端。

**API：**

```text
def find_nearest_transcript(gtf: str | Path, contig: str, position: int, *, end: int | None=None, strand: str | None=None, cache_dir: str | Path | None=None, verbose: bool=True) -> dict[str, Any]:
```

**输入：** 签名列出必选参数、可选参数与原生默认值；参数含义和输出如下。

```text
Find the transcript whose locus is nearest to a position or interval, even if none overlaps
(``Genome.nearest_transcript``).

``distance`` is 0 when the transcript overlaps the position (``position`` to ``end``).
```

**输出：** 返回可 JSON 序列化的字典；文件输出返回路径。可选依赖缺失时抛出带安装提示的 BackendUnavailableError。

## `dnakit.annotation.pyensembl_p36_functions.get_annotation_by_id`

**作用：** Get the annotation of one gene, transcript or exon by ID (``Genome.gene_by_id``, ``transcript_by_id``,

**处理方法：** 由模块中的显式接口校验输入并调用对应后端。

**API：**

```text
def get_annotation_by_id(gtf: str | Path, feature: str, identifier: str, *, cache_dir: str | Path | None=None, verbose: bool=True) -> dict[str, Any]:
```

**输入：** 签名列出必选参数、可选参数与原生默认值；参数含义和输出如下。

```text
Get the annotation of one gene, transcript or exon by ID (``Genome.gene_by_id``, ``transcript_by_id``,
``exon_by_id``; the name and locus look-ups ``gene_name_of_gene_id``, ``locus_of_*_id`` and
``transcript_name_of_transcript_id`` are included in the record).

``feature`` is ``gene``, ``transcript`` or ``exon``.
```

**输出：** 返回可 JSON 序列化的字典；文件输出返回路径。可选依赖缺失时抛出带安装提示的 BackendUnavailableError。

## `dnakit.annotation.pyensembl_p36_functions.list_annotation_contigs`

**作用：** List the contig names that carry genes in an annotation (``Genome.contigs``).

**处理方法：** 由模块中的显式接口校验输入并调用对应后端。

**API：**

```text
def list_annotation_contigs(gtf: str | Path, *, cache_dir: str | Path | None=None, verbose: bool=True) -> dict[str, Any]:
```

**输入：** 签名列出必选参数、可选参数与原生默认值；参数含义和输出如下。

```text
List the contig names that carry genes in an annotation (``Genome.contigs``).
```

**输出：** 返回可 JSON 序列化的字典；文件输出返回路径。可选依赖缺失时抛出带安装提示的 BackendUnavailableError。

## `dnakit.annotation.pyensembl_p36_functions.list_annotation_features`

**作用：** List genes, transcripts, exons or protein IDs, filtered by contig, strand or biotype

**处理方法：** 由模块中的显式接口校验输入并调用对应后端。

**API：**

```text
def list_annotation_features(gtf: str | Path, feature: str='gene', *, contig: str | None=None, strand: str | None=None, biotype: str | None=None, output: str='ids', cache_dir: str | Path | None=None, max_results: int=1000, verbose: bool=True) -> dict[str, Any]:
```

**输入：** 签名列出必选参数、可选参数与原生默认值；参数含义和输出如下。

```text
List genes, transcripts, exons or protein IDs, filtered by contig, strand or biotype
(``Genome.gene_ids`` / ``gene_names`` / ``genes``, ``transcript_ids`` / ``transcript_names`` /
``transcripts``, ``exon_ids`` / ``exons``, ``protein_ids``).

``feature`` is ``gene``, ``transcript``, ``exon`` or ``protein``; ``output`` is ``ids``, ``names``
(genes and transcripts) or ``records`` (genes, transcripts, exons). ``biotype`` applies to genes and
transcripts only (pyensembl has no exon or protein biotype).
```

**输出：** 返回可 JSON 序列化的字典；文件输出返回路径。可选依赖缺失时抛出带安装提示的 BackendUnavailableError。

## `dnakit.annotation.pyensembl_p36_functions.related_annotation_features`

**作用：** Follow the relations between genes, transcripts, exons and proteins in an annotation.

**处理方法：** 由模块中的显式接口校验输入并调用对应后端。

**API：**

```text
def related_annotation_features(gtf: str | Path, from_kind: str, identifier: str, to_kind: str, *, cache_dir: str | Path | None=None, verbose: bool=True) -> dict[str, Any]:
```

**输入：** 签名列出必选参数、可选参数与原生默认值；参数含义和输出如下。

```text
Follow the relations between genes, transcripts, exons and proteins in an annotation.

``from_kind`` is ``gene_id``, ``gene_name``, ``transcript_id``, ``transcript_name``, ``exon_id`` or
``protein_id``; ``to_kind`` is ``gene``, ``transcript`` or ``exon``. Supported pairs: protein to gene
(``gene_by_protein_id``) and to transcript (``transcript_by_protein_id``); transcript to gene
(``gene_name_of_transcript_id`` / ``_name``) and to exons (``exon_ids_of_transcript_id`` / ``_name``);
transcript name to transcripts (``transcripts_by_name``); exon to gene (``gene_name_of_exon_id``) and to
transcripts (``transcript_ids_of_exon_id``); gene to transcripts (``transcript_ids_of_gene_id`` /
``_gene_name``) and to exons (``exon_ids_of_gene_id`` / ``_gene_name``).
```

**输出：** 返回可 JSON 序列化的字典；文件输出返回路径。可选依赖缺失时抛出带安装提示的 BackendUnavailableError。

## `dnakit.annotation.pyensembl_p36_functions.map_position_to_transcript_offset`

**作用：** Convert a genomic position to the 0-based offset in the spliced transcript sequence

**处理方法：** 由模块中的显式接口校验输入并调用对应后端。

**API：**

```text
def map_position_to_transcript_offset(gtf: str | Path, transcript_id: str, position: int, *, cache_dir: str | Path | None=None, verbose: bool=True) -> dict[str, Any]:
```

**输入：** 签名列出必选参数、可选参数与原生默认值；参数含义和输出如下。

```text
Convert a genomic position to the 0-based offset in the spliced transcript sequence
(``Transcript.spliced_offset``).

The offset counts exon bases only, from the transcript's 5' end (strand-aware). A position that lies
in an intron or outside every exon is reported with ``offset: null`` and the reason.
```

**输出：** 返回可 JSON 序列化的字典；文件输出返回路径。可选依赖缺失时抛出带安装提示的 BackendUnavailableError。

## `dnakit.annotation.pyensembl_p36_functions.genomic_locus_offset`

**作用：** Offset of a position (or a position range) from the stranded start of a locus

**处理方法：** 由模块中的显式接口校验输入并调用对应后端。

**API：**

```text
def genomic_locus_offset(contig: str, start: int, end: int, strand: str, *, position: int | None=None, range_start: int | None=None, range_end: int | None=None, verbose: bool=True) -> dict[str, Any]:
```

**输入：** 签名列出必选参数、可选参数与原生默认值；参数含义和输出如下。

```text
Offset of a position (or a position range) from the stranded start of a locus
(``Locus.offset``, ``Locus.offset_range``).

On the plus strand the offset of position 13 in 10..20 is 3, on the minus strand 7. Give ``position``
or both ``range_start`` and ``range_end``; positions must lie inside the locus.
```

**输出：** 返回可 JSON 序列化的字典；文件输出返回路径。可选依赖缺失时抛出带安装提示的 BackendUnavailableError。

## `dnakit.annotation.pyensembl_p36_functions.compare_genomic_loci`

**作用：** Distance and overlap length of two loci (``Locus.distance_to_locus`` / ``distance_to_interval`` and

**处理方法：** 由模块中的显式接口校验输入并调用对应后端。

**API：**

```text
def compare_genomic_loci(first: Mapping[str, Any], second: Mapping[str, Any], *, ignore_strand: bool=False, verbose: bool=True) -> dict[str, Any]:
```

**输入：** 签名列出必选参数、可选参数与原生默认值；参数含义和输出如下。

```text
Distance and overlap length of two loci (``Locus.distance_to_locus`` / ``distance_to_interval`` and
``Locus.overlap_length``).

Each locus is ``{"contig", "start", "end", "strand"}`` (1-based inclusive). ``distance`` is 0 for
overlapping loci and ``null`` for loci on different contigs (or strands, unless ``ignore_strand``);
``overlap_length`` is the number of shared bases (0 when they cannot overlap).
```

**输出：** 返回可 JSON 序列化的字典；文件输出返回路径。可选依赖缺失时抛出带安装提示的 BackendUnavailableError。

## 可执行示例

**示例代码：**

```python
from dnakit.annotation.pyensembl_p36_functions import genomic_locus_offset
r = genomic_locus_offset('chr1', 10, 20, '-', position=13, verbose=False)
print(r['offset'])
```

**示例输出：**

```text
7
```
