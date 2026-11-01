# Chamaeleo DNA 存储编码

本页记录 P37 相关的真实 DNAKit 接口。函数按所需可选依赖延迟导入；坐标、默认值与参数选择以各接口说明为准。

验证范围说明：Chamaeleo 1.34 was tested from a wheel extracted into scratch; it is absent from the default validation environment.

## `dnakit.molecular_lab.chamaeleo_p37_functions.chamaeleo_encode_direct`

**作用：** Encode bit segments into DNA with the 2-bit direct mapping (00 to A, 01 to C, 10 to G, 11 to T).

**处理方法：** 由模块中的显式接口校验输入并调用对应后端。

**API：**

```text
def chamaeleo_encode_direct(bit_segments: Sequence[Sequence[int] | str], *, verbose: bool=True) -> dict[str, Any]:
```

**输入：** 签名列出必选参数、可选参数与原生默认值；参数含义和输出如下。

```text
Encode bit segments into DNA with the 2-bit direct mapping (00 to A, 01 to C, 10 to G, 11 to T).

Segment lengths must be even. The result holds the DNA strings and the information density in bits per nucleotide.
```

**输出：** 返回可 JSON 序列化的字典；文件输出返回路径。可选依赖缺失时抛出带安装提示的 BackendUnavailableError。

## `dnakit.molecular_lab.chamaeleo_p37_functions.chamaeleo_decode_direct`

**作用：** Decode DNA sequences back into bit segments with the 2-bit direct mapping.

**处理方法：** 由模块中的显式接口校验输入并调用对应后端。

**API：**

```text
def chamaeleo_decode_direct(dna_sequences: Sequence[str], *, verbose: bool=True) -> dict[str, Any]:
```

**输入：** 签名列出必选参数、可选参数与原生默认值；参数含义和输出如下。

```text
Decode DNA sequences back into bit segments with the 2-bit direct mapping.
```

**输出：** 返回可 JSON 序列化的字典；文件输出返回路径。可选依赖缺失时抛出带安装提示的 BackendUnavailableError。

## `dnakit.molecular_lab.chamaeleo_p37_functions.chamaeleo_encode_church`

**作用：** Encode bits with Church et al. 2012: bit 0 is A or C, bit 1 is G or T, and a fourth identical base in a row is avoided.

**处理方法：** 由模块中的显式接口校验输入并调用对应后端。

**API：**

```text
def chamaeleo_encode_church(bit_segments: Sequence[Sequence[int] | str], *, seed: int | None=None, verbose: bool=True) -> dict[str, Any]:
```

**输入：** 签名列出必选参数、可选参数与原生默认值；参数含义和输出如下。

```text
Encode bits with Church et al. 2012: bit 0 is A or C, bit 1 is G or T, and a fourth identical base in a row is avoided.

The choice between the two bases is random; pass ``seed`` for a reproducible result.
```

**输出：** 返回可 JSON 序列化的字典；文件输出返回路径。可选依赖缺失时抛出带安装提示的 BackendUnavailableError。

## `dnakit.molecular_lab.chamaeleo_p37_functions.chamaeleo_decode_church`

**作用：** Decode Church et al. 2012 DNA (A or C is 0, G or T is 1) back into bit segments.

**处理方法：** 由模块中的显式接口校验输入并调用对应后端。

**API：**

```text
def chamaeleo_decode_church(dna_sequences: Sequence[str], *, verbose: bool=True) -> dict[str, Any]:
```

**输入：** 签名列出必选参数、可选参数与原生默认值；参数含义和输出如下。

```text
Decode Church et al. 2012 DNA (A or C is 0, G or T is 1) back into bit segments.
```

**输出：** 返回可 JSON 序列化的字典；文件输出返回路径。可选依赖缺失时抛出带安装提示的 BackendUnavailableError。

## `dnakit.molecular_lab.chamaeleo_p37_functions.chamaeleo_encode_goldman`

**作用：** Encode bytes with Goldman et al. 2013: a ternary Huffman code whose trits pick the next base by rotation.

**处理方法：** 由模块中的显式接口校验输入并调用对应后端。

**API：**

```text
def chamaeleo_encode_goldman(bit_segments: Sequence[Sequence[int] | str], *, fixed_huffman: bool=True, support_nucleotide: str='A', verbose: bool=True) -> dict[str, Any]:
```

**输入：** 签名列出必选参数、可选参数与原生默认值；参数含义和输出如下。

```text
Encode bytes with Goldman et al. 2013: a ternary Huffman code whose trits pick the next base by rotation.

Segment lengths must be multiples of 8. ``fixed_huffman=True`` uses the published code. ``False`` builds an adaptive
tree from the data (``chamaeleo_goldman_huffman_tree``), uses it for encoding and returns it as ``huffman_tree``;
``chamaeleo_decode_goldman`` needs it. (Chamaeleo's own adaptive path passes the argument ``3`` as the data size
and fails for most inputs, so the tree is built here with the intended size and multiple.)
```

**输出：** 返回可 JSON 序列化的字典；文件输出返回路径。可选依赖缺失时抛出带安装提示的 BackendUnavailableError。

## `dnakit.molecular_lab.chamaeleo_p37_functions.chamaeleo_decode_goldman`

**作用：** Decode Goldman et al. 2013 DNA into bit segments; pass the ``huffman_tree`` returned by an adaptive encoding (default: the published code).

**处理方法：** 由模块中的显式接口校验输入并调用对应后端。

**API：**

```text
def chamaeleo_decode_goldman(dna_sequences: Sequence[str], *, support_nucleotide: str='A', huffman_tree: Sequence[str | None] | None=None, segment_length: int | None=None, verbose: bool=True) -> dict[str, Any]:
```

**输入：** 签名列出必选参数、可选参数与原生默认值；参数含义和输出如下。

```text
Decode Goldman et al. 2013 DNA into bit segments; pass the ``huffman_tree`` returned by an adaptive encoding (default: the published code).
```

**输出：** 返回可 JSON 序列化的字典；文件输出返回路径。可选依赖缺失时抛出带安装提示的 BackendUnavailableError。

## `dnakit.molecular_lab.chamaeleo_p37_functions.chamaeleo_goldman_huffman_tree`

**作用：** Build the adaptive ternary Huffman code of Goldman coding from the byte frequencies of the data (256 codes, ``None`` for unused bytes).

**处理方法：** 由模块中的显式接口校验输入并调用对应后端。

**API：**

```text
def chamaeleo_goldman_huffman_tree(bit_segments: Sequence[Sequence[int] | str], *, multiple: int=3, verbose: bool=True) -> dict[str, Any]:
```

**输入：** 签名列出必选参数、可选参数与原生默认值；参数含义和输出如下。

```text
Build the adaptive ternary Huffman code of Goldman coding from the byte frequencies of the data (256 codes, ``None`` for unused bytes).
```

**输出：** 返回可 JSON 序列化的字典；文件输出返回路径。可选依赖缺失时抛出带安装提示的 BackendUnavailableError。

## `dnakit.molecular_lab.chamaeleo_p37_functions.chamaeleo_encode_grass`

**作用：** Encode 16-bit blocks into three base-47 codewords (Grass et al. 2015, 9 nucleotides per block).

**处理方法：** 由模块中的显式接口校验输入并调用对应后端。

**API：**

```text
def chamaeleo_encode_grass(bit_segments: Sequence[Sequence[int] | str], *, base_values: Sequence[int] | None=None, verbose: bool=True) -> dict[str, Any]:
```

**输入：** 签名列出必选参数、可选参数与原生默认值；参数含义和输出如下。

```text
Encode 16-bit blocks into three base-47 codewords (Grass et al. 2015, 9 nucleotides per block).

Segment lengths must be multiples of 16. ``base_values`` is an optional custom list of 48 base values.
```

**输出：** 返回可 JSON 序列化的字典；文件输出返回路径。可选依赖缺失时抛出带安装提示的 BackendUnavailableError。

## `dnakit.molecular_lab.chamaeleo_p37_functions.chamaeleo_decode_grass`

**作用：** Decode Grass et al. 2015 DNA (9 nucleotides per 16 bits) into bit segments.

**处理方法：** 由模块中的显式接口校验输入并调用对应后端。

**API：**

```text
def chamaeleo_decode_grass(dna_sequences: Sequence[str], *, base_values: Sequence[int] | None=None, segment_length: int | None=None, verbose: bool=True) -> dict[str, Any]:
```

**输入：** 签名列出必选参数、可选参数与原生默认值；参数含义和输出如下。

```text
Decode Grass et al. 2015 DNA (9 nucleotides per 16 bits) into bit segments.
```

**输出：** 返回可 JSON 序列化的字典；文件输出返回路径。可选依赖缺失时抛出带安装提示的 BackendUnavailableError。

## `dnakit.molecular_lab.chamaeleo_p37_functions.chamaeleo_encode_blawat`

**作用：** Encode 8-bit blocks into 5 nucleotides (Blawat et al. 2016) under GC-content and homopolymer constraints.

**处理方法：** 由模块中的显式接口校验输入并调用对应后端。

**API：**

```text
def chamaeleo_encode_blawat(bit_segments: Sequence[Sequence[int] | str], *, seed: int | None=None, verbose: bool=True) -> dict[str, Any]:
```

**输入：** 签名列出必选参数、可选参数与原生默认值；参数含义和输出如下。

```text
Encode 8-bit blocks into 5 nucleotides (Blawat et al. 2016) under GC-content and homopolymer constraints.

Segment lengths must be multiples of 8. Where several codewords satisfy the constraints one is chosen at random; set ``seed``.
```

**输出：** 返回可 JSON 序列化的字典；文件输出返回路径。可选依赖缺失时抛出带安装提示的 BackendUnavailableError。

## `dnakit.molecular_lab.chamaeleo_p37_functions.chamaeleo_decode_blawat`

**作用：** Decode Blawat et al. 2016 DNA (5 nucleotides per byte) into bit segments.

**处理方法：** 由模块中的显式接口校验输入并调用对应后端。

**API：**

```text
def chamaeleo_decode_blawat(dna_sequences: Sequence[str], *, segment_length: int | None=None, verbose: bool=True) -> dict[str, Any]:
```

**输入：** 签名列出必选参数、可选参数与原生默认值；参数含义和输出如下。

```text
Decode Blawat et al. 2016 DNA (5 nucleotides per byte) into bit segments.
```

**输出：** 返回可 JSON 序列化的字典；文件输出返回路径。可选依赖缺失时抛出带安装提示的 BackendUnavailableError。

## `dnakit.molecular_lab.chamaeleo_p37_functions.chamaeleo_encode_dna_fountain`

**作用：** Encode bit segments as DNA Fountain droplets (Erlich and Zielinski 2017) screened by homopolymer and GC limits.

**处理方法：** 由模块中的显式接口校验输入并调用对应后端。

**API：**

```text
def chamaeleo_encode_dna_fountain(bit_segments: Sequence[Sequence[int] | str], *, homopolymer: int=4, gc_bias: float=0.2, redundancy: float=0.07, header_size: int=4, c_dist: float=0.1, delta: float=0.05, need_pre_check: bool=False, seed: int | None=None, verbose: bool=True) -> dict[str, Any]:
```

**输入：** 签名列出必选参数、可选参数与原生默认值；参数含义和输出如下。

```text
Encode bit segments as DNA Fountain droplets (Erlich and Zielinski 2017) screened by homopolymer and GC limits.

Droplets are generated until ``ceil(segments * (1 + redundancy))`` valid sequences exist. Each starts with a
``header_size``-byte seed. Segment lengths must be even (odd ones get a leading 0). ``need_pre_check=True`` test-decodes first.
Returns ``decode_packets``, the segment count that ``chamaeleo_decode_dna_fountain`` needs.
```

**输出：** 返回可 JSON 序列化的字典；文件输出返回路径。可选依赖缺失时抛出带安装提示的 BackendUnavailableError。

## `dnakit.molecular_lab.chamaeleo_p37_functions.chamaeleo_decode_dna_fountain`

**作用：** Decode DNA Fountain droplets by peeling; ``decode_packets`` is the original segment count and must match the encoding parameters.

**处理方法：** 由模块中的显式接口校验输入并调用对应后端。

**API：**

```text
def chamaeleo_decode_dna_fountain(dna_sequences: Sequence[str], decode_packets: int, *, header_size: int=4, c_dist: float=0.1, delta: float=0.05, segment_length: int | None=None, verbose: bool=True) -> dict[str, Any]:
```

**输入：** 签名列出必选参数、可选参数与原生默认值；参数含义和输出如下。

```text
Decode DNA Fountain droplets by peeling; ``decode_packets`` is the original segment count and must match the encoding parameters.
```

**输出：** 返回可 JSON 序列化的字典；文件输出返回路径。可选依赖缺失时抛出带安装提示的 BackendUnavailableError。

## `dnakit.molecular_lab.chamaeleo_p37_functions.chamaeleo_encode_yin_yang`

**作用：** Encode bit segments with the YinYang code, which merges two bit segments into one DNA sequence.

**处理方法：** 由模块中的显式接口校验输入并调用对应后端。

**API：**

```text
def chamaeleo_encode_yin_yang(bit_segments: Sequence[Sequence[int] | str], *, faster: bool=False, max_iterations: int=100, max_ratio: float=0.8, max_homopolymer: int=4, max_content: float=0.6, virtual_nucleotide: str='A', seed: int | None=None, verbose: bool=True) -> dict[str, Any]:
```

**输入：** 签名列出必选参数、可选参数与原生默认值；参数含义和输出如下。

```text
Encode bit segments with the YinYang code, which merges two bit segments into one DNA sequence.

``faster=False`` is the normal setting (good and bad segments are paired deliberately); ``faster=True`` pairs segments at
random and may add more filler sequences. Sequences obey ``max_homopolymer`` and ``max_content`` (GC limit).
Returns ``total_count``, the segment count that ``chamaeleo_decode_yin_yang`` needs to drop filler sequences.
```

**输出：** 返回可 JSON 序列化的字典；文件输出返回路径。可选依赖缺失时抛出带安装提示的 BackendUnavailableError。

## `dnakit.molecular_lab.chamaeleo_p37_functions.chamaeleo_decode_yin_yang`

**作用：** Decode YinYang DNA into bit segments; filler segments with an index above ``total_count`` are dropped.

**处理方法：** 由模块中的显式接口校验输入并调用对应后端。

**API：**

```text
def chamaeleo_decode_yin_yang(dna_sequences: Sequence[str], total_count: int, *, virtual_nucleotide: str='A', segment_length: int | None=None, verbose: bool=True) -> dict[str, Any]:
```

**输入：** 签名列出必选参数、可选参数与原生默认值；参数含义和输出如下。

```text
Decode YinYang DNA into bit segments; filler segments with an index above ``total_count`` are dropped.

The decoded segments come in pairs (upper then lower bit segment) and are not sorted; the index prefix of each segment
(the first ``bit_length(total_count)`` bits) gives its position when the data was encoded with an index.
```

**输出：** 返回可 JSON 序列化的字典；文件输出返回路径。可选依赖缺失时抛出带安装提示的 BackendUnavailableError。

## `dnakit.molecular_lab.chamaeleo_p37_functions.chamaeleo_screen_dna`

**作用：** Screen one DNA sequence for homopolymer runs and GC content (``screen.homopolymer``, ``screen.gc_content``, ``screen.check``).

**处理方法：** 由模块中的显式接口校验输入并调用对应后端。

**API：**

```text
def chamaeleo_screen_dna(sequence: str, *, max_homopolymer: int | None=None, max_gc_deviation: float | None=None, verbose: bool=True) -> dict[str, Any]:
```

**输入：** 签名列出必选参数、可选参数与原生默认值；参数含义和输出如下。

```text
Screen one DNA sequence for homopolymer runs and GC content (``screen.homopolymer``, ``screen.gc_content``, ``screen.check``).

The homopolymer test fails when a base repeats more than ``max_homopolymer`` times in a row. The GC test uses
Chamaeleo's ``max_content`` ``m`` (here ``0.5 + max_gc_deviation``) and passes when ``1 - m <= GC fraction <= m``.
Give either limit or both (the combined ``check``); a limit left as None is not tested.
```

**输出：** 返回可 JSON 序列化的字典；文件输出返回路径。可选依赖缺失时抛出带安装提示的 BackendUnavailableError。

## `dnakit.molecular_lab.chamaeleo_p37_functions.chamaeleo_read_dna_file`

**作用：** Read DNA sequences from a text file with one sequence per line (``read_dna_file``).

**处理方法：** 由模块中的显式接口校验输入并调用对应后端。

**API：**

```text
def chamaeleo_read_dna_file(path: str | Path, *, max_sequences: int=1000000, verbose: bool=True) -> dict[str, Any]:
```

**输入：** 签名列出必选参数、可选参数与原生默认值；参数含义和输出如下。

```text
Read DNA sequences from a text file with one sequence per line (``read_dna_file``).
```

**输出：** 返回可 JSON 序列化的字典；文件输出返回路径。可选依赖缺失时抛出带安装提示的 BackendUnavailableError。

## `dnakit.molecular_lab.chamaeleo_p37_functions.chamaeleo_write_dna_file`

**作用：** Write DNA sequences to a text file, one per line (``write_dna_file``). The file must not exist yet.

**处理方法：** 由模块中的显式接口校验输入并调用对应后端。

**API：**

```text
def chamaeleo_write_dna_file(sequences: Sequence[str], output_path: str | Path, *, verbose: bool=True) -> dict[str, Any]:
```

**输入：** 签名列出必选参数、可选参数与原生默认值；参数含义和输出如下。

```text
Write DNA sequences to a text file, one per line (``write_dna_file``). The file must not exist yet.
```

**输出：** 返回可 JSON 序列化的字典；文件输出返回路径。可选依赖缺失时抛出带安装提示的 BackendUnavailableError。

## `dnakit.molecular_lab.chamaeleo_p37_functions.chamaeleo_transcode_file_to_dna`

**作用：** Transcode a digital file to DNA (``TranscodePipeline.transcode(direction='t_c')``).

**处理方法：** 由模块中的显式接口校验输入并调用对应后端。

**API：**

```text
def chamaeleo_transcode_file_to_dna(input_path_: str | Path, output_path: str | Path, *, scheme: str='direct', scheme_options: Mapping[str, Any] | None=None, error_correction: str | None=None, check_size: int=3, segment_length: int=120, index: bool=True, index_length: int | None=None, seed: int | None=None, verbose: bool=True) -> dict[str, Any]:
```

**输入：** 签名列出必选参数、可选参数与原生默认值；参数含义和输出如下。

```text
Transcode a digital file to DNA (``TranscodePipeline.transcode(direction='t_c')``).

The file bytes are cut into ``segment_length``-bit segments (the last one zero-padded), numbered with an ``index``
prefix, protected with ``error_correction`` (``hamming`` or ``reed_solomon`` with ``check_size`` bytes), encoded with
``scheme`` (one of direct, church, goldman, grass, blawat, dna_fountain, yin_yang; ``scheme_options`` are its constructor
options) and written one sequence per line. The result holds ``decode_parameters`` for ``chamaeleo_transcode_dna_to_file``.
```

**输出：** 返回可 JSON 序列化的字典；文件输出返回路径。可选依赖缺失时抛出带安装提示的 BackendUnavailableError。

## `dnakit.molecular_lab.chamaeleo_p37_functions.chamaeleo_transcode_dna_to_file`

**作用：** Transcode DNA (one sequence per line) back to a digital file (``TranscodePipeline.transcode(direction='t_s')``).

**处理方法：** 由模块中的显式接口校验输入并调用对应后端。

**API：**

```text
def chamaeleo_transcode_dna_to_file(input_path_: str | Path, output_path: str | Path, *, bit_size: int, segment_length: int, scheme: str='direct', scheme_options: Mapping[str, Any] | None=None, error_correction: str | None=None, check_size: int=3, index: bool=True, index_length: int | None=None, payload_count: int | None=None, verbose: bool=True) -> dict[str, Any]:
```

**输入：** 签名列出必选参数、可选参数与原生默认值；参数含义和输出如下。

```text
Transcode DNA (one sequence per line) back to a digital file (``TranscodePipeline.transcode(direction='t_s')``).

``bit_size`` (the original file size in bits) and the other parameters are the ``decode_parameters`` returned by
``chamaeleo_transcode_file_to_dna``; ``segment_length`` is the length of an encoded bit segment including the index
and error-correction bits (not the payload length). ``dna_fountain`` also needs ``decode_packets`` and ``yin_yang`` needs ``total_count`` in
``scheme_options`` (``payload_count`` fills these in when they are missing).
```

**输出：** 返回可 JSON 序列化的字典；文件输出返回路径。可选依赖缺失时抛出带安装提示的 BackendUnavailableError。

## `dnakit.molecular_lab.chamaeleo_p37_functions.chamaeleo_evaluate_robustness`

**作用：** Simulate DNA-storage errors and measure how much of each file is recovered (``RobustnessPipeline.evaluate``).

**处理方法：** 由模块中的显式接口校验输入并调用对应后端。

**API：**

```text
def chamaeleo_evaluate_robustness(file_paths: Mapping[str, str | Path], coding_schemes: Sequence[str], *, error_corrections: Sequence[str | None]=(None,), check_size: int=3, needed_indices: Sequence[bool] | None=None, scheme_options: Mapping[str, Mapping[str, Any]] | None=None, nucleotide_insertion: float=0.0, nucleotide_mutation: float=0.0, nucleotide_deletion: float=0.0, sequence_loss: float=0.0, iterations: int=1, segment_length: int=120, index_length: int | None=None, seed: int | None=None, log_path: str | Path | None=None, verbose: bool=True) -> dict[str, Any]:
```

**输入：** 签名列出必选参数、可选参数与原生默认值；参数含义和输出如下。

```text
Simulate DNA-storage errors and measure how much of each file is recovered (``RobustnessPipeline.evaluate``).

Every file is encoded with every scheme and error-correction option, a fraction of the sequences is lost, and the
given rates of nucleotide insertions, mutations and deletions (0 to 1, per sequence count) are applied before decoding.
``needed_indices`` (one flag per scheme, default all True) says whether a scheme gets an index. Each task reports the
information density, runtimes and the recovery success rate per iteration. ``log_path`` (a new file) receives the CSV log
that ``chamaeleo_choose_best_scheme`` reads.
```

**输出：** 返回可 JSON 序列化的字典；文件输出返回路径。可选依赖缺失时抛出带安装提示的 BackendUnavailableError。

## `dnakit.molecular_lab.chamaeleo_p37_functions.chamaeleo_profile_coding_features`

**作用：** Count the GC-content and homopolymer-length distributions of the DNA each coding scheme makes from each file (``BasicFeaturePipeline.calculate``).

**处理方法：** 由模块中的显式接口校验输入并调用对应后端。

**API：**

```text
def chamaeleo_profile_coding_features(file_paths: Mapping[str, str | Path], coding_schemes: Sequence[str], *, needed_indices: Sequence[bool] | None=None, scheme_options: Mapping[str, Mapping[str, Any]] | None=None, segment_length: int=120, index_length: int | None=None, seed: int | None=None, log_path: str | Path | None=None, verbose: bool=True) -> dict[str, Any]:
```

**输入：** 签名列出必选参数、可选参数与原生默认值；参数含义和输出如下。

```text
Count the GC-content and homopolymer-length distributions of the DNA each coding scheme makes from each file (``BasicFeaturePipeline.calculate``).

``gc_distribution`` has 101 counts (sequences per rounded GC percent) and ``homopolymer_distribution`` counts
sequences by their longest homopolymer run. ``log_path`` (a new file) receives the log for ``chamaeleo_choose_best_scheme``.
```

**输出：** 返回可 JSON 序列化的字典；文件输出返回路径。可选依赖缺失时抛出带安装提示的 BackendUnavailableError。

## `dnakit.molecular_lab.chamaeleo_p37_functions.chamaeleo_choose_best_scheme`

**作用：** Compare coding schemes from the logs of ``chamaeleo_profile_coding_features`` and ``chamaeleo_evaluate_robustness`` (``OptimalChoicePipeline.calculate_best``).

**处理方法：** 由模块中的显式接口校验输入并调用对应后端。

**API：**

```text
def chamaeleo_choose_best_scheme(feature_log_path: str | Path, robustness_log_path: str | Path, *, verbose: bool=True) -> dict[str, Any]:
```

**输入：** 签名列出必选参数、可选参数与原生默认值；参数含义和输出如下。

```text
Compare coding schemes from the logs of ``chamaeleo_profile_coding_features`` and ``chamaeleo_evaluate_robustness`` (``OptimalChoicePipeline.calculate_best``).

Returns one row per scheme with the average information density, the share of sequences with 40 to 60 percent GC, the
maximum GC bias, the maximum homopolymer length and the average recovery rate, as Chamaeleo's table prints them.
Chamaeleo's parser reads logs written with the default 101-bin GC distribution.
```

**输出：** 返回可 JSON 序列化的字典；文件输出返回路径。可选依赖缺失时抛出带安装提示的 BackendUnavailableError。

## 可执行示例

**示例代码：**

```python
from dnakit.molecular_lab.chamaeleo_p37_functions import chamaeleo_encode_direct, chamaeleo_decode_direct
r = chamaeleo_encode_direct([[0, 0, 0, 1, 1, 0, 1, 1]], verbose=False)
print(r['dna'])
print(chamaeleo_decode_direct(r['dna'], verbose=False)['bit_strings'])
```

**示例输出：**

```text
['ACGT']
['00011011']
```
