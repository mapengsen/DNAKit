# HGVS 校验、映射与数据提供方

本页记录 P35 相关的真实 DNAKit 接口。函数按所需可选依赖延迟导入；坐标、默认值与参数选择以各接口说明为准。

验证范围说明：Mapping supports local reference snapshots and caller-selected UTA access. NCBI PostgreSQL mirror methods are excluded because the required database is not publicly provided.

## `dnakit.variants.hgvs_p35_functions.parse_hgvs_posedit`

**作用：** Parse the position-plus-edit part of an HGVS string, for example ``22+1A>T`` (``c``) or ``Ala22Trp`` (``p``).

**处理方法：** 由模块中的显式接口校验输入并调用对应后端。

**API：**

```text
def parse_hgvs_posedit(text: str, *, variant_type: str='c', verbose: bool=True) -> dict[str, Any]:
```

**输入：** 签名列出必选参数、可选参数与原生默认值；参数含义和输出如下。

```text
Parse the position-plus-edit part of an HGVS string, for example ``22+1A>T`` (``c``) or ``Ala22Trp`` (``p``).

``variant_type`` selects the grammar rule: ``c`` (coding DNA, ``Parser.parse_c_posedit``), ``g``
(genomic), ``m`` (mitochondrial), ``n`` (non-coding), ``p`` (protein) or ``r`` (RNA). The result gives
the normalized text, the position (start, end) and the edit (type, reference, alternate).
```

**输出：** 返回可 JSON 序列化的字典；文件输出返回路径。可选依赖缺失时抛出带安装提示的 BackendUnavailableError。

## `dnakit.variants.hgvs_p35_functions.parse_hgvs_position`

**作用：** Parse an HGVS position or interval string, for example ``NM_000001.1:c.26+2_57-3`` or ``26+2_57-3``.

**处理方法：** 由模块中的显式接口校验输入并调用对应后端。

**API：**

```text
def parse_hgvs_position(text: str, *, variant_type: str | None=None, verbose: bool=True) -> dict[str, Any]:
```

**输入：** 签名列出必选参数、可选参数与原生默认值；参数含义和输出如下。

```text
Parse an HGVS position or interval string, for example ``NM_000001.1:c.26+2_57-3`` or ``26+2_57-3``.

A string with an accession and type prefix (``NM_000001.1:c.26+2_57-3``) is parsed by
``parse_hgvs_position``, or by ``parse_c_hgvs_position`` and its siblings when ``variant_type`` (``c``,
``g``, ``m``, ``n``, ``p`` or ``r``) is given. A bare interval (``26+2_57-3``) needs ``variant_type`` and
uses ``parse_c_interval`` and its siblings. The result gives the start and end positions, with base,
offset and datum for coding positions.
```

**输出：** 返回可 JSON 序列化的字典；文件输出返回路径。可选依赖缺失时抛出带安装提示的 BackendUnavailableError。

## `dnakit.variants.hgvs_p35_functions.validate_hgvs_position`

**作用：** Validate an HGVS position or interval intrinsically (no reference data needed).

**处理方法：** 由模块中的显式接口校验输入并调用对应后端。

**API：**

```text
def validate_hgvs_position(text: str, *, variant_type: str='c', verbose: bool=True) -> dict[str, Any]:
```

**输入：** 签名列出必选参数、可选参数与原生默认值；参数含义和输出如下。

```text
Validate an HGVS position or interval intrinsically (no reference data needed).

``text`` is an interval such as ``26+2_57-3`` of the given ``variant_type``. The check is hgvs'
``Interval.validate``: base zero, bases below one for sequence-start positions and start after end are
errors. Returns ``valid`` with the validation level and message.
```

**输出：** 返回可 JSON 序列化的字典；文件输出返回路径。可选依赖缺失时抛出带安装提示的 BackendUnavailableError。

## `dnakit.variants.hgvs_p35_functions.validate_hgvs_variant`

**作用：** Validate an HGVS variant: intrinsic, extrinsic or both.

**处理方法：** 由模块中的显式接口校验输入并调用对应后端。

**API：**

```text
def validate_hgvs_variant(variant: str, *, level: str='intrinsic', strict: bool=True, data_source: str='local', local_data: dict[str, Any] | None=None, uta_db_url: str | None=None, verbose: bool=True) -> dict[str, Any]:
```

**输入：** 签名列出必选参数、可选参数与原生默认值；参数含义和输出如下。

```text
Validate an HGVS variant: intrinsic, extrinsic or both.

``level="intrinsic"`` checks the string itself (accession and type agree, insertion and deletion
lengths, position order; ``IntrinsicValidator``) and needs no data. ``"extrinsic"`` checks the variant
against the reference sequence and transcript bounds (``ExtrinsicValidator``; needs ``data_source``).
``"full"`` does both (``Validator``, as ``hgvs.easy.validate``). ``strict`` treats warnings as failures,
as hgvs does by default. Returns ``valid`` and, if invalid, the hgvs message.
```

**输出：** 返回可 JSON 序列化的字典；文件输出返回路径。可选依赖缺失时抛出带安装提示的 BackendUnavailableError。

## `dnakit.variants.hgvs_p35_functions.fill_hgvs_reference`

**作用：** Fill the reference bases of a deletion, duplication, delins or repeat variant from the reference sequence.

**处理方法：** 由模块中的显式接口校验输入并调用对应后端。

**API：**

```text
def fill_hgvs_reference(variant: str, *, data_source: str='local', local_data: dict[str, Any] | None=None, uta_db_url: str | None=None, genomic_accession: str | None=None, alt_aln_method: str='splign', verbose: bool=True) -> dict[str, Any]:
```

**输入：** 签名列出必选参数、可选参数与原生默认值；参数含义和输出如下。

```text
Fill the reference bases of a deletion, duplication, delins or repeat variant from the reference sequence.

``NM_x:c.5_7del`` becomes ``NM_x:c.5_7delACG`` (``SequenceVariant.fill_ref``; the text is written with the
reference bases shown, and ``reference`` holds them). For transcript variants
the reference is read through the transcript's genomic alignment (``genomic_accession`` selects it);
variants that carry their reference already are returned unchanged.
```

**输出：** 返回可 JSON 序列化的字典；文件输出返回路径。可选依赖缺失时抛出带安装提示的 BackendUnavailableError。

## `dnakit.variants.hgvs_p35_functions.hgvs_length_change`

**作用：** Get the net length change of an HGVS variant (``PosEdit.length_change``).

**处理方法：** 由模块中的显式接口校验输入并调用对应后端。

**API：**

```text
def hgvs_length_change(variant: str, *, verbose: bool=True) -> dict[str, Any]:
```

**输入：** 签名列出必选参数、可选参数与原生默认值；参数含义和输出如下。

```text
Get the net length change of an HGVS variant (``PosEdit.length_change``).

Negative values are net deletions, zero is a balanced change such as an SNV, positive values are net
insertions. Deletion lengths can be inferred from exact interval coordinates even
when deleted bases are omitted. Undefined changes raise ``ConfigurationError``.
```

**输出：** 返回可 JSON 序列化的字典；文件输出返回路径。可选依赖缺失时抛出带安装提示的 BackendUnavailableError。

## `dnakit.variants.hgvs_p35_functions.map_hgvs_variant`

**作用：** Map an HGVS variant between genomic (g.), coding (c.) and non-coding (n.) coordinates.

**处理方法：** 由模块中的显式接口校验输入并调用对应后端。

**API：**

```text
def map_hgvs_variant(variant: str, *, to_type: str, transcript_accession: str | None=None, genomic_accession: str | None=None, assembly_name: str | None=None, alt_aln_method: str='splign', prevalidation_level: str | None='EXTRINSIC', normalize: bool=True, replace_reference: bool=True, data_source: str='local', local_data: dict[str, Any] | None=None, uta_db_url: str | None=None, verbose: bool=True) -> dict[str, Any]:
```

**输入：** 签名列出必选参数、可选参数与原生默认值；参数含义和输出如下。

```text
Map an HGVS variant between genomic (g.), coding (c.) and non-coding (n.) coordinates.

One function for the hgvs mappers ``g_to_c``, ``g_to_n``, ``g_to_t``, ``n_to_g``, ``c_to_g``,
``t_to_g``, ``c_to_n`` and ``n_to_c``: the input type is read from the variant and ``to_type`` is ``g``,
``c``, ``n`` or ``t`` (``t`` picks c. or n. by transcript type). To a transcript type from g. give
``transcript_accession``; to g. from a transcript give ``genomic_accession`` or ``assembly_name`` (then
hgvs' ``AssemblyMapper`` finds the chromosome and normalizes, as ``hgvs.easy`` does); c. <-> n. needs
neither. ``prevalidation_level`` (``EXTRINSIC`` by default, as in hgvs), ``normalize`` and
``replace_reference`` follow hgvs' defaults. ``data_source`` is ``local`` (``local_data``) or ``uta``.
```

**输出：** 返回可 JSON 序列化的字典；文件输出返回路径。可选依赖缺失时抛出带安装提示的 BackendUnavailableError。

## `dnakit.variants.hgvs_p35_functions.map_hgvs_interval`

**作用：** Map a position interval between genomic, coding and non-coding coordinates with an ``AlignmentMapper``.

**处理方法：** 由模块中的显式接口校验输入并调用对应后端。

**API：**

```text
def map_hgvs_interval(interval: str, *, from_type: str, to_type: str, transcript_accession: str, genomic_accession: str | None=None, alt_aln_method: str='splign', strict_bounds: bool | None=None, data_source: str='local', local_data: dict[str, Any] | None=None, uta_db_url: str | None=None, verbose: bool=True) -> dict[str, Any]:
```

**输入：** 签名列出必选参数、可选参数与原生默认值；参数含义和输出如下。

```text
Map a position interval between genomic, coding and non-coding coordinates with an ``AlignmentMapper``.

``interval`` is a bare interval such as ``26+2_57-3`` of ``from_type``; ``from_type`` and ``to_type``
are two of ``g``, ``c`` and ``n``. The mapper is built for ``transcript_accession`` and its
``genomic_accession`` (not needed for c. <-> n.). ``strict_bounds`` overrides hgvs' default of requiring
transcript positions inside the transcript. This is the interval form of ``map_hgvs_variant``
(``TranscriptMapper``/``AlignmentMapper`` ``g_to_c`` and siblings).
```

**输出：** 返回可 JSON 序列化的字典；文件输出返回路径。可选依赖缺失时抛出带安装提示的 BackendUnavailableError。

## `dnakit.variants.hgvs_p35_functions.predict_hgvs_protein_consequence`

**作用：** Predict the protein consequence (p.) of a coding (c.) variant.

**处理方法：** 由模块中的显式接口校验输入并调用对应后端。

**API：**

```text
def predict_hgvs_protein_consequence(variant: str, *, protein_accession: str | None=None, genomic_accession: str | None=None, assembly_name: str | None=None, alt_aln_method: str='splign', data_source: str='local', local_data: dict[str, Any] | None=None, uta_db_url: str | None=None, verbose: bool=True) -> dict[str, Any]:
```

**输入：** 签名列出必选参数、可选参数与原生默认值；参数含义和输出如下。

```text
Predict the protein consequence (p.) of a coding (c.) variant.

This is ``VariantMapper.c_to_p`` (``AssemblyMapper.c_to_p`` when ``assembly_name`` is given and
``AssemblyMapper.t_to_p``-style handling: a non-coding n. variant returns ``non-coding``). The protein
accession defaults to the transcript's own; hgvs uses the standard translation table, adds
``p.(...)`` uncertainty brackets for inferred changes (``inferred_p_is_uncertain``) and reports
``p.?``, ``p.=`` or frameshift/extension notation as appropriate.
```

**输出：** 返回可 JSON 序列化的字典；文件输出返回路径。可选依赖缺失时抛出带安装提示的 BackendUnavailableError。

## `dnakit.variants.hgvs_p35_functions.project_hgvs_variant`

**作用：** Project a coding variant or interval between two transcripts through a common genomic sequence.

**处理方法：** 由模块中的显式接口校验输入并调用对应后端。

**API：**

```text
def project_hgvs_variant(variant: str, *, source_accession: str, destination_accession: str, genomic_accession: str, direction: str='forward', source_alt_aln_method: str='splign', destination_alt_aln_method: str='splign', data_source: str='local', local_data: dict[str, Any] | None=None, uta_db_url: str | None=None, verbose: bool=True) -> dict[str, Any]:
```

**输入：** 签名列出必选参数、可选参数与原生默认值；参数含义和输出如下。

```text
Project a coding variant or interval between two transcripts through a common genomic sequence.

``direction="forward"`` projects from ``source_accession`` to ``destination_accession``
(``Projector.project_variant_forward`` / ``project_interval_forward``); ``"backward"`` projects a variant
on the destination back to the source. ``variant`` is a full c. variant on the transcript of the
starting side, or a bare c. interval such as ``26+2_57-3``. The transcripts must be on the same strand
of ``genomic_accession``.
```

**输出：** 返回可 JSON 序列化的字典；文件输出返回路径。可选依赖缺失时抛出带安装提示的 BackendUnavailableError。

## `dnakit.variants.hgvs_p35_functions.list_hgvs_relevant_transcripts`

**作用：** List the transcripts whose genomic alignment covers a genomic (g.) variant.

**处理方法：** 由模块中的显式接口校验输入并调用对应后端。

**API：**

```text
def list_hgvs_relevant_transcripts(variant: str, *, assembly_name: str='GRCh38', alt_aln_method: str='splign', data_source: str='local', local_data: dict[str, Any] | None=None, uta_db_url: str | None=None, verbose: bool=True) -> dict[str, Any]:
```

**输入：** 签名列出必选参数、可选参数与原生默认值；参数含义和输出如下。

```text
List the transcripts whose genomic alignment covers a genomic (g.) variant.

This is ``AssemblyMapper.relevant_transcripts`` (UTA's ``get_tx_for_region`` query): the transcripts
whose aligned span starts before the variant start and reaches its end, on the variant's chromosome.
```

**输出：** 返回可 JSON 序列化的字典；文件输出返回路径。可选依赖缺失时抛出带安装提示的 BackendUnavailableError。

## `dnakit.variants.hgvs_p35_functions.hgvs_to_vcf_allele`

**作用：** Convert a genomic (g.) HGVS variant to a VCF-style allele (chromosome, position, ref, alt).

**处理方法：** 由模块中的显式接口校验输入并调用对应后端。

**API：**

```text
def hgvs_to_vcf_allele(variant: str, *, assembly_name: str='GRCh38', data_source: str='local', local_data: dict[str, Any] | None=None, uta_db_url: str | None=None, verbose: bool=True) -> dict[str, Any]:
```

**输入：** 签名列出必选参数、可选参数与原生默认值；参数含义和输出如下。

```text
Convert a genomic (g.) HGVS variant to a VCF-style allele (chromosome, position, ref, alt).

This is hgvs' experimental ``Babelfish.hgvs_to_vcf``: the variant is left-normalized, deletions and
insertions are anchored on the preceding base, duplications become an insertion, and the chromosome
uses the official name (no ``chr`` prefix) of ``assembly_name``. Reference bases are read through the
data source.
```

**输出：** 返回可 JSON 序列化的字典；文件输出返回路径。可选依赖缺失时抛出带安装提示的 BackendUnavailableError。

## `dnakit.variants.hgvs_p35_functions.fetch_hgvs_reference_sequence`

**作用：** Fetch a reference sequence by accession and zero-based interval ``[start_i, end_i)``.

**处理方法：** 由模块中的显式接口校验输入并调用对应后端。

**API：**

```text
def fetch_hgvs_reference_sequence(accession: str, *, start_i: int | None=None, end_i: int | None=None, data_source: str='local', local_data: dict[str, Any] | None=None, uta_db_url: str | None=None, verbose: bool=True) -> dict[str, Any]:
```

**输入：** 签名列出必选参数、可选参数与原生默认值；参数含义和输出如下。

```text
Fetch a reference sequence by accession and zero-based interval ``[start_i, end_i)``.

``data_source="local"`` reads ``local_data['sequences']``; ``"uta"`` goes through the UTA data
provider's sequence fetcher, which reads SeqRepo (``HGVS_SEQREPO_DIR`` or ``HGVS_SEQREPO_URL``) or
falls back to NCBI/Ensembl over the network. Intervals are as in ``SeqFetcher.fetch_seq``.
```

**输出：** 返回可 JSON 序列化的字典；文件输出返回路径。可选依赖缺失时抛出带安装提示的 BackendUnavailableError。

## `dnakit.variants.hgvs_p35_functions.query_hgvs_transcript_database`

**作用：** Query the transcript database behind hgvs (UTA schema) for genes, transcripts and alignments.

**处理方法：** 由模块中的显式接口校验输入并调用对应后端。

**API：**

```text
def query_hgvs_transcript_database(query: str, *, gene: str | None=None, transcript_accession: str | None=None, genomic_accession: str | None=None, alt_aln_method: str | None=None, assembly_name: str | None=None, data_source: str='local', local_data: dict[str, Any] | None=None, uta_db_url: str | None=None, verbose: bool=True) -> dict[str, Any]:
```

**输入：** 签名列出必选参数、可选参数与原生默认值；参数含义和输出如下。

```text
Query the transcript database behind hgvs (UTA schema) for genes, transcripts and alignments.

``query`` selects the data-provider call: ``gene_info`` (HGNC symbol ``gene``; ``get_gene_info``),
``tx_exons`` (``transcript_accession``, ``genomic_accession``, ``alt_aln_method``; exon structure on the
genome alignment; ``get_tx_exons``), ``tx_for_gene`` (transcripts of ``gene``; ``get_tx_for_gene``),
``tx_identity_info`` (transcript lengths and CDS bounds; ``get_tx_identity_info``), ``tx_info``
(alignment info for transcript, genomic accession and method; ``get_tx_info``), ``tx_mapping_options``
(alignment options of ``transcript_accession``; ``get_tx_mapping_options``), ``similar_transcripts``
(transcripts with identical CDS or exon structure; ``get_similar_transcripts``), ``protein_accession``
(the protein of a transcript; ``get_pro_ac_for_tx_ac``) and ``assembly_map`` (the chromosome accessions
of ``assembly_name`` and their names; ``get_assembly_map``, which reads hgvs' bundled assembly data).
Rows are returned as dictionaries.
```

**输出：** 返回可 JSON 序列化的字典；文件输出返回路径。可选依赖缺失时抛出带安装提示的 BackendUnavailableError。

## 可执行示例

**示例代码：**

```python
from dnakit.variants.hgvs_p35_functions import hgvs_length_change, parse_hgvs_position
print(hgvs_length_change('NM_1.1:c.5_7del', verbose=False)['length_change'])
print(parse_hgvs_position('26+2_57-3', variant_type='c', verbose=False)['start']['datum'])
```

**示例输出：**

```text
-3
CDS_START
```
