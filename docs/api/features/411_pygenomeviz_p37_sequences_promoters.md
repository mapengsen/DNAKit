# pyGenomeViz 序列、比对与启动子

本页记录 P37 相关的真实 DNAKit 接口。函数按所需可选依赖延迟导入；坐标、默认值与参数选择以各接口说明为准。

验证范围说明：Promoters use the native 1.x API or equivalent matplotlib geometry on 0.4 feature tracks. Genbank.write_cds_fasta exports proteins upstream; the inventory label claiming nucleotide output is incorrect.

## `dnakit.visualization.pygenomeviz_p37_functions.read_contig_sequences`

**作用：** Read the contig (record) sequences of a GenBank or FASTA file as an ID to sequence listing.

**处理方法：** 由模块中的显式接口校验输入并调用对应后端。

**API：**

```text
def read_contig_sequences(source: str | Path, *, file_format: str='auto', include_sequences: bool=True, verbose: bool=True) -> dict[str, Any]:
```

**输入：** 签名列出必选参数、可选参数与原生默认值；参数含义和输出如下。

```text
Read the contig (record) sequences of a GenBank or FASTA file as an ID to sequence listing.

Covers ``Genbank.get_seqid2seq`` (``file_format='genbank'``) and ``Fasta.get_seqid2seq`` /
``get_seqid2record`` (``'fasta'``). FASTA files with duplicate IDs are rejected, as pyGenomeViz does.
```

**输出：** 返回可 JSON 序列化的字典；文件输出返回路径。可选依赖缺失时抛出带安装提示的 BackendUnavailableError。

## `dnakit.visualization.pygenomeviz_p37_functions.read_gff_contig_sizes`

**作用：** Contig lengths from GFF ``##sequence-region`` headers (``Gff.get_seqid2size``).

**处理方法：** 由模块中的显式接口校验输入并调用对应后端。

**API：**

```text
def read_gff_contig_sizes(source: str | Path, *, verbose: bool=True) -> dict[str, Any]:
```

**输入：** 签名列出必选参数、可选参数与原生默认值；参数含义和输出如下。

```text
Contig lengths from GFF ``##sequence-region`` headers (``Gff.get_seqid2size``).

A contig without a header gets the largest end coordinate of its features, which can be shorter than the contig.
```

**输出：** 返回可 JSON 序列化的字典；文件输出返回路径。可选依赖缺失时抛出带安装提示的 BackendUnavailableError。

## `dnakit.visualization.pygenomeviz_p37_functions.list_contig_features`

**作用：** List the features of a type and strand for every contig of a GenBank or GFF file.

**处理方法：** 由模块中的显式接口校验输入并调用对应后端。

**API：**

```text
def list_contig_features(source: str | Path, *, file_format: str='auto', feature_type: str | Sequence[str] | None='CDS', strand: int | None=None, verbose: bool=True) -> dict[str, Any]:
```

**输入：** 签名列出必选参数、可选参数与原生默认值；参数含义和输出如下。

```text
List the features of a type and strand for every contig of a GenBank or GFF file.

Covers ``Genbank.get_seqid2features`` and ``Gff.get_seqid2features``. ``feature_type=None`` keeps every type.
GenBank features that wrap over the genome origin are dropped. ``strand`` is 1, -1 or None (any).
Coordinates are 0-based half-open.
```

**输出：** 返回可 JSON 序列化的字典；文件输出返回路径。可选依赖缺失时抛出带安装提示的 BackendUnavailableError。

## `dnakit.visualization.pygenomeviz_p37_functions.extract_genome_features`

**作用：** Extract features of one contig, optionally inside a range (``Genbank.extract_features``, ``Gff.extract_features``).

**处理方法：** 由模块中的显式接口校验输入并调用对应后端。

**API：**

```text
def extract_genome_features(source: str | Path, *, file_format: str='auto', feature_type: str | Sequence[str] | None='CDS', strand: int | None=None, target_range: Sequence[int] | None=None, seqid: str | None=None, verbose: bool=True) -> dict[str, Any]:
```

**输入：** 签名列出必选参数、可选参数与原生默认值；参数含义和输出如下。

```text
Extract features of one contig, optionally inside a range (``Genbank.extract_features``, ``Gff.extract_features``).

GenBank uses the first record that has a matching feature and GFF the first seqid, unless ``seqid`` is given.
``target_range`` is ``(start, end)``: GenBank compares 0-based half-open feature coordinates, GFF the 1-based
record coordinates, as pyGenomeViz does.
```

**输出：** 返回可 JSON 序列化的字典；文件输出返回路径。可选依赖缺失时抛出带安装提示的 BackendUnavailableError。

## `dnakit.visualization.pygenomeviz_p37_functions.extract_gff_exon_features`

**作用：** Build exon-structured features from the ``ID``/``Parent`` relation of a feature type and its exons (``Gff.extract_exon_features``).

**处理方法：** 由模块中的显式接口校验输入并调用对应后端。

**API：**

```text
def extract_gff_exon_features(source: str | Path, *, feature_type: str='mRNA', strand: int | None=None, target_range: Sequence[int] | None=None, seqid: str | None=None, verbose: bool=True) -> dict[str, Any]:
```

**输入：** 签名列出必选参数、可选参数与原生默认值；参数含义和输出如下。

```text
Build exon-structured features from the ``ID``/``Parent`` relation of a feature type and its exons (``Gff.extract_exon_features``).

A parent with one exon keeps its own span; with two or more exons it gets a ``parts`` list of exon
``[start, end]`` pairs (0-based half-open, sorted); a parent without exons is skipped. Exons are matched to the
most recent parent of ``feature_type`` in file order, as pyGenomeViz does.
```

**输出：** 返回可 JSON 序列化的字典；文件输出返回路径。可选依赖缺失时抛出带安装提示的 BackendUnavailableError。

## `dnakit.visualization.pygenomeviz_p37_functions.select_features_within_range`

**作用：** Keep the features that lie completely inside ``(start, end)`` (``pygenomeviz.utils.extract_features_within_range``).

**处理方法：** 由模块中的显式接口校验输入并调用对应后端。

**API：**

```text
def select_features_within_range(features: Sequence[Mapping[str, Any]], target_range: Sequence[int], *, verbose: bool=True) -> dict[str, Any]:
```

**输入：** 签名列出必选参数、可选参数与原生默认值；参数含义和输出如下。

```text
Keep the features that lie completely inside ``(start, end)`` (``pygenomeviz.utils.extract_features_within_range``).

Each feature is a mapping with integer ``start`` and ``end`` (for example rows from ``extract_genome_features``).
```

**输出：** 返回可 JSON 序列化的字典；文件输出返回路径。可选依赖缺失时抛出带安装提示的 BackendUnavailableError。

## `dnakit.visualization.pygenomeviz_p37_functions.write_genome_fasta`

**作用：** Write the genome sequences of a GenBank or FASTA file as one-line-per-sequence FASTA (``write_genome_fasta``).

**处理方法：** 由模块中的显式接口校验输入并调用对应后端。

**API：**

```text
def write_genome_fasta(source: str | Path, output_path: str | Path, *, file_format: str='auto', verbose: bool=True) -> dict[str, Any]:
```

**输入：** 签名列出必选参数、可选参数与原生默认值；参数含义和输出如下。

```text
Write the genome sequences of a GenBank or FASTA file as one-line-per-sequence FASTA (``write_genome_fasta``).
```

**输出：** 返回可 JSON 序列化的字典；文件输出返回路径。可选依赖缺失时抛出带安装提示的 BackendUnavailableError。

## `dnakit.visualization.pygenomeviz_p37_functions.write_genbank_cds_fasta`

**作用：** Write the CDS features of a GenBank file with a translation to protein FASTA (``Genbank.write_cds_fasta``).

**处理方法：** 由模块中的显式接口校验输入并调用对应后端。

**API：**

```text
def write_genbank_cds_fasta(source: str | Path, output_path: str | Path, *, verbose: bool=True) -> dict[str, Any]:
```

**输入：** 签名列出必选参数、可选参数与原生默认值；参数含义和输出如下。

```text
Write the CDS features of a GenBank file with a translation to protein FASTA (``Genbank.write_cds_fasta``).

Each record is named ``GENE<n>_<protein_id>|<seqid>|<start>_<end>_<strand>|`` with the product as description.
CDS features without a ``/translation`` and features that wrap the origin are skipped. For the nucleotide
sequence of every feature use ``extract_gffutils_feature_sequence`` or Biopython.
```

**输出：** 返回可 JSON 序列化的字典；文件输出返回路径。可选依赖缺失时抛出带安装提示的 BackendUnavailableError。

## `dnakit.visualization.pygenomeviz_p37_functions.calc_genbank_gc_content`

**作用：** GC content (percent) in a sliding window along a GenBank genome (``Genbank.calc_gc_content``).

**处理方法：** 由模块中的显式接口校验输入并调用对应后端。

**API：**

```text
def calc_genbank_gc_content(source: str | Path, *, window_size: int | None=None, step_size: int | None=None, seqid: str | None=None, verbose: bool=True) -> dict[str, Any]:
```

**输入：** 签名列出必选参数、可选参数与原生默认值；参数含义和输出如下。

```text
GC content (percent) in a sliding window along a GenBank genome (``Genbank.calc_gc_content``).

Positions run from 0 to the sequence length in steps of ``step_size`` (default length/1000); each value is the GC
percentage of the window centred on the position (default window length/500). The first record is used unless ``seqid`` is given.
```

**输出：** 返回可 JSON 序列化的字典；文件输出返回路径。可选依赖缺失时抛出带安装提示的 BackendUnavailableError。

## `dnakit.visualization.pygenomeviz_p37_functions.parse_alignment_coordinates`

**作用：** Parse a genome alignment result file into alignment coordinates (``AlignCoord.parse_*_file``).

**处理方法：** 由模块中的显式接口校验输入并调用对应后端。

**API：**

```text
def parse_alignment_coordinates(alignment_file: str | Path, *, file_format: str='blast', query_id: str='query', ref_id: str='ref', seqtype: str='nucleotide', names: Sequence[str] | None=None, ref_index: int=0, verbose: bool=True) -> dict[str, Any]:
```

**输入：** 签名列出必选参数、可选参数与原生默认值；参数含义和输出如下。

```text
Parse a genome alignment result file into alignment coordinates (``AlignCoord.parse_*_file``).

``file_format`` is ``blast`` (BLAST or LAST tabular, ``-outfmt 6``), ``mummer`` (``show-coords -H -T -r -k`` of nucmer,
``seqtype='nucleotide'``, or promer, ``'protein'``) or ``progressive_mauve`` (backbone file; ``names`` lists the genomes,
blocks shorter than 20 bp or missing in two genomes are dropped). Completely overlapping blocks are removed.
Coordinates are converted to 0-based starts.
```

**输出：** 返回可 JSON 序列化的字典；文件输出返回路径。可选依赖缺失时抛出带安装提示的 BackendUnavailableError。

## `dnakit.visualization.pygenomeviz_p37_functions.write_alignment_coordinates_tsv`

**作用：** Write alignment coordinates as the 12-column pyGenomeViz TSV (``AlignCoord.write``); unknown scores are ``na``.

**处理方法：** 由模块中的显式接口校验输入并调用对应后端。

**API：**

```text
def write_alignment_coordinates_tsv(coords: Sequence[Mapping[str, Any]], output_path: str | Path, *, verbose: bool=True) -> dict[str, Any]:
```

**输入：** 签名列出必选参数、可选参数与原生默认值；参数含义和输出如下。

```text
Write alignment coordinates as the 12-column pyGenomeViz TSV (``AlignCoord.write``); unknown scores are ``na``.
```

**输出：** 返回可 JSON 序列化的字典；文件输出返回路径。可选依赖缺失时抛出带安装提示的 BackendUnavailableError。

## `dnakit.visualization.pygenomeviz_p37_functions.read_alignment_coordinates_tsv`

**作用：** Read alignment coordinates from the pyGenomeViz TSV written by ``write_alignment_coordinates_tsv`` (``AlignCoord.read``).

**处理方法：** 由模块中的显式接口校验输入并调用对应后端。

**API：**

```text
def read_alignment_coordinates_tsv(coords_file: str | Path, *, verbose: bool=True) -> dict[str, Any]:
```

**输入：** 签名列出必选参数、可选参数与原生默认值；参数含义和输出如下。

```text
Read alignment coordinates from the pyGenomeViz TSV written by ``write_alignment_coordinates_tsv`` (``AlignCoord.read``).
```

**输出：** 返回可 JSON 序列化的字典；文件输出返回路径。可选依赖缺失时抛出带安装提示的 BackendUnavailableError。

## `dnakit.visualization.pygenomeviz_p37_functions.filter_alignment_coordinates`

**作用：** Filter alignment coordinates by length, identity and e-value (``AlignCoord.filter``) and optionally drop contained blocks.

**处理方法：** 由模块中的显式接口校验输入并调用对应后端。

**API：**

```text
def filter_alignment_coordinates(coords: Sequence[Mapping[str, Any]], *, length_thr: int | None=None, identity_thr: float | None=None, evalue_thr: float | None=None, remove_overlapping: bool=False, verbose: bool=True) -> dict[str, Any]:
```

**输入：** 签名列出必选参数、可选参数与原生默认值；参数含义和输出如下。

```text
Filter alignment coordinates by length, identity and e-value (``AlignCoord.filter``) and optionally drop contained blocks.

A block is dropped when either of its lengths is below ``length_thr``, its identity is below ``identity_thr`` or its
e-value is above ``evalue_thr``; a threshold of ``None`` or 0 disables the test, and blocks without a score pass the
score tests. ``remove_overlapping=True`` also removes blocks completely inside another block of the same
query/reference pair (``AlignCoord.filter_overlap``).
```

**输出：** 返回可 JSON 序列化的字典；文件输出返回路径。可选依赖缺失时抛出带安装提示的 BackendUnavailableError。

## `dnakit.visualization.pygenomeviz_p37_functions.align_genomes_blast`

**作用：** Align each genome FASTA with the next one by BLAST (``Blast.run``, ``pgv-blast``).

**处理方法：** 由模块中的显式接口校验输入并调用对应后端。

**API：**

```text
def align_genomes_blast(seq_files: Sequence[str | Path], output_dir: str | Path, *, seqtype: str='nucleotide', evalue: float=0.001, threads: int=1, timeout_seconds: float=3600.0, blast_dir: str | Path | None=None, verbose: bool=True) -> dict[str, Any]:
```

**输入：** 签名列出必选参数、可选参数与原生默认值；参数含义和输出如下。

```text
Align each genome FASTA with the next one by BLAST (``Blast.run``, ``pgv-blast``).

``seqtype='nucleotide'`` runs ``blastn`` and ``'protein'`` runs ``tblastx`` (translated comparison). For every
neighbouring pair the second file becomes a ``makeblastdb`` database, the first is the query, and the tabular
result (``-outfmt 6 -evalue <evalue>``) is parsed into alignment coordinates with contained hits removed.
```

**输出：** 返回可 JSON 序列化的字典；文件输出返回路径。可选依赖缺失时抛出带安装提示的 BackendUnavailableError。

## `dnakit.visualization.pygenomeviz_p37_functions.align_genomes_last`

**作用：** Align each genome FASTA with the next one by LAST (``Last.run``), following the LAST cookbook genome recipe.

**处理方法：** 由模块中的显式接口校验输入并调用对应后端。

**API：**

```text
def align_genomes_last(seq_files: Sequence[str | Path], output_dir: str | Path, *, threads: int=1, timeout_seconds: float=3600.0, last_dir: str | Path | None=None, verbose: bool=True) -> dict[str, Any]:
```

**输入：** 签名列出必选参数、可选参数与原生默认值；参数含义和输出如下。

```text
Align each genome FASTA with the next one by LAST (``Last.run``), following the LAST cookbook genome recipe.

For every neighbouring pair: ``lastdb`` on the second file, ``lastal -D 1e9 --split-f=MAF+`` for the first, ``last-split -r``
to make the alignment one-to-one, and ``maf-convert -n blasttab``; the tabular result is parsed into alignment coordinates.
```

**输出：** 返回可 JSON 序列化的字典；文件输出返回路径。可选依赖缺失时抛出带安装提示的 BackendUnavailableError。

## `dnakit.visualization.pygenomeviz_p37_functions.draw_genome_segment_features`

**作用：** Draw gene glyphs on one genome segment with pyGenomeViz (``add_feature(s)``, ``add_exon_feature(s)``) and save the figure.

**处理方法：** 由模块中的显式接口校验输入并调用对应后端。

**API：**

```text
def draw_genome_segment_features(features: Sequence[Mapping[str, Any]], output_path: str | Path, *, segment_name: str='genome', segment_size: int, fig_width: float=15.0, dpi: int=150, verbose: bool=True) -> dict[str, Any]:
```

**输入：** 签名列出必选参数、可选参数与原生默认值；参数含义和输出如下。

```text
Draw gene glyphs on one genome segment with pyGenomeViz (``add_feature(s)``, ``add_exon_feature(s)``) and save the figure.

Each feature is a mapping. A gene arrow needs ``start``, ``end`` and ``strand`` (1 or -1); a gene with exons
gives ``parts`` (a list of ``[start, end]``) instead of ``start``/``end``, and is drawn as exon boxes joined by
intron lines. Optional keys are ``label``, ``facecolor``, ``edgecolor``, ``plotstyle`` (``bigarrow``, ``arrow``,
``box``, ``bigbox``, ``rbox``, ``bigrbox``) and ``linewidth``. Needs pyGenomeViz 0.4.x and matplotlib.
The promoter glyph of pyGenomeViz 1.x is not available.
```

**输出：** 返回可 JSON 序列化的字典；文件输出返回路径。可选依赖缺失时抛出带安装提示的 BackendUnavailableError。

## `dnakit.visualization.pygenomeviz_p37_functions.export_genome_html_viewer`

**作用：** Export GenBank genomes (CDS features as tracks) as one interactive HTML viewer with hover tooltips (``GenomeViz.savefig_html``).

**处理方法：** 由模块中的显式接口校验输入并调用对应后端。

**API：**

```text
def export_genome_html_viewer(gbk_resources: Sequence[str | Path], html_path: str | Path, *, fig_width: float=15.0, verbose: bool=True) -> dict[str, Any]:
```

**输入：** 签名列出必选参数、可选参数与原生默认值；参数含义和输出如下。

```text
Export GenBank genomes (CDS features as tracks) as one interactive HTML viewer with hover tooltips (``GenomeViz.savefig_html``).

Needs pyGenomeViz 0.4.x. The output is a single HTML file that embeds the SVG figure and the viewer script.
```

**输出：** 返回可 JSON 序列化的字典；文件输出返回路径。可选依赖缺失时抛出带安装提示的 BackendUnavailableError。

## `dnakit.visualization.pygenomeviz_p37_promoter_functions.draw_genome_promoters`

**作用：** Draw promoter arrows using FeatureSegment.add_promoter and return image path plus glyph data.

**处理方法：** 由模块中的显式接口校验输入并调用对应后端。

**API：**

```text
def draw_genome_promoters(promoters: Sequence[Mapping[str, Any]], output_path: str, *, segment_size: int, segment_name: str='genome', verbose: bool=True) -> dict[str, Any]:
```

**输入：** 签名列出必选参数、可选参数与原生默认值；参数含义和输出如下。

```text
Draw promoter arrows using FeatureSegment.add_promoter and return image path plus glyph data.

Each promoter has ``position`` (0-based), optional ``strand`` (+1/-1),
``length`` (relative to segment, default 0.02), and ``color`` (default black).
Uses FeatureSegment.add_promoter on pyGenomeViz 1.x. On 0.4 the same
right-angle promoter geometry is drawn on the feature track's matplotlib
axes, preserving compatibility with existing DNAKit genome plots.
```

**输出：** 返回可 JSON 序列化的字典；文件输出返回路径。可选依赖缺失时抛出带安装提示的 BackendUnavailableError。

## 可执行示例

**示例代码：**

```python
from pathlib import Path
from dnakit.visualization.pygenomeviz_p37_functions import read_contig_sequences
from dnakit.visualization.pygenomeviz_p37_promoter_functions import draw_genome_promoters
Path('seq.fa').write_text('>chr1\nACGT\n')
print(read_contig_sequences('seq.fa', verbose=False)['contigs'][0]['sequence'])
r = draw_genome_promoters([{'position': 100, 'strand': 1}], 'promoter.svg', segment_size=1000, verbose=False)
print(r['promoters'][0]['arrow_end'], Path(r['image_path']).is_file())
```

**示例输出：**

```text
ACGT
[120.0, 1.0] True
```
