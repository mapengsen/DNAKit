# DnaChisel 设计规则与 DendroPy 群体模型

本页记录 P38 相关的真实 DNAKit 接口。函数按所需可选依赖延迟导入；坐标、默认值与参数选择以各接口说明为准。

验证范围说明：BLAST/Bowtie/primer3-dependent rules need their corresponding programs or packages and caller databases. No external program or remote resource is invoked by the examples.

## `dnakit.molecular_lab.dnachisel_p38_functions.evaluate_dnachisel_design_rule`

**作用：** Score a DNA sequence against one DnaChisel design rule, and optionally optimise it for the rule.

**处理方法：** 由模块中的显式接口校验输入并调用对应后端。

**API：**

```text
def evaluate_dnachisel_design_rule(sequence: str, rule: str, *, parameters: Mapping[str, Any] | None=None, circular: bool=False, role: str='objective', optimize: bool=False, verbose: bool=True) -> dict[str, Any]:
```

**输入：** 签名列出必选参数、可选参数与原生默认值；参数含义和输出如下。

```text
Score a DNA sequence against one DnaChisel design rule, and optionally optimise it for the rule.

``rule`` selects the specification: ``AllowPrimer`` (primer-friendly region: Tm window, no repeats, homologies
or heterodimers), ``AvoidBlastMatches`` (BLAST+ ``blast_db`` or ``sequences``), ``AvoidHeterodimerization``
(``other_primers_sequences``, ``tmax``), ``AvoidMatches`` (``match_length`` and a Bowtie ``bowtie_index``),
``EnforceChoice`` (``choices``), ``EnforcePatternOccurence`` (``pattern``, ``occurences``),
``EnforceRegionsCompatibility`` (``locations`` and ``condition`` ``different`` or ``min_hamming_distance`` with
``min_distance``), ``EnforceSequence`` (degenerate ``sequence`` at ``location``), ``EnforceTerminalGCContent``
(``window_size``, ``mini``, ``maxi``), ``SequenceLengthBounds`` (``min_length``, ``max_length``),
``UniquifyAllKmers`` (``k``), ``AvoidRareCodons`` (``min_frequency``, ``species``), ``MatchTargetCodonUsage``
(``species``) and ``HarmonizeRCA`` (``original_species``, ``species``). ``parameters`` holds the DnaChisel
keyword arguments of that rule; ``location`` is ``[start, end]`` or ``[start, end, strand]``.

``role`` is ``objective`` (the sequence is only scored; default) or ``constraint`` (DnaChisel first resolves
constraints that restrict nucleotides, which can change the sequence). With ``optimize=True`` the rule is used
as an objective (or constraint) and the sequence is optimised with DnaChisel's random-mutation optimiser; the
result holds the evaluation before and after and the new sequence. Optimisation is stochastic in DnaChisel.
```

**输出：** 返回可 JSON 序列化的字典；文件输出返回路径。可选依赖缺失时抛出带安装提示的 BackendUnavailableError。

## `dnakit.molecular_lab.dnachisel_p38_functions.build_pssm_pattern_from_sequences`

**作用：** Build a PSSM pattern from same-length DNA sequences with pseudocounts, using

**处理方法：** 由模块中的显式接口校验输入并调用对应后端。

**API：**

```text
def build_pssm_pattern_from_sequences(sequences: Sequence[str], *, name: str='unnamed', pseudocounts: Any='jaspar', threshold: float | None=None, relative_threshold: float | None=None, scan_sequence: str | None=None, verbose: bool=True) -> dict[str, Any]:
```

**输入：** 签名列出必选参数、可选参数与原生默认值；参数含义和输出如下。

```text
Build a PSSM pattern from same-length DNA sequences with pseudocounts, using
``dnachisel.MotifPssmPattern.from_sequences`` (and ``apply_pseudocounts``).

``pseudocounts`` is ``jaspar`` (automatic, Biopython), a mapping ``{"A": .., "C": .., "G": .., "T": ..}`` or
null for none. A match threshold is given absolute (``threshold``) or relative (``relative_threshold``, 0 to
1 between the minimum and maximum possible score). With ``scan_sequence`` the forward-strand matches are
returned. Returns the log-odds matrix by base.
```

**输出：** 返回可 JSON 序列化的字典；文件输出返回路径。可选依赖缺失时抛出带安装提示的 BackendUnavailableError。

## `dnakit.molecular_lab.dnachisel_p38_functions.read_pssm_patterns_from_motif_file`

**作用：** Read motifs from a JASPAR, MEME or TRANSFAC file and return them as PSSM patterns, using

**处理方法：** 由模块中的显式接口校验输入并调用对应后端。

**API：**

```text
def read_pssm_patterns_from_motif_file(motifs_path: str, file_format: str, *, pseudocounts: Any='jaspar', threshold: float | None=None, relative_threshold: float | None=None, scan_sequence: str | None=None, verbose: bool=True) -> dict[str, Any]:
```

**输入：** 签名列出必选参数、可选参数与原生默认值；参数含义和输出如下。

```text
Read motifs from a JASPAR, MEME or TRANSFAC file and return them as PSSM patterns, using
``dnachisel.MotifPssmPattern.list_from_file`` (Bio.motifs parsers).

``file_format`` is ``jaspar``, ``meme`` (MEME XML, as read by Biopython) or ``TRANSFAC``. Pseudocounts and thresholds are as in
``build_pssm_pattern_from_sequences``; with ``scan_sequence`` the forward-strand matches of every motif are
listed.
```

**输出：** 返回可 JSON 序列化的字典；文件输出返回路径。可选依赖缺失时抛出带安装提示的 BackendUnavailableError。

## `dnakit.phylogeny.dendropy_p38_functions.read_multi_phylip_alignments`

**作用：** Read a multi-block PHYLIP file (several alignments in one file) with dendropy's ``multiphylip`` reader.

**处理方法：** 由模块中的显式接口校验输入并调用对应后端。

**API：**

```text
def read_multi_phylip_alignments(path: str, *, data_type: str='dna', strict: bool=False, interleaved: bool=False, limit: int=100, verbose: bool=True) -> dict[str, Any]:
```

**输入：** 签名列出必选参数、可选参数与原生默认值；参数含义和输出如下。

```text
Read a multi-block PHYLIP file (several alignments in one file) with dendropy's ``multiphylip`` reader.

Every block starts with its own ``ntax nchar`` header line. ``data_type`` is ``dna`` or ``rna``; ``strict``
selects strict (10-character names) rather than relaxed PHYLIP; ``interleaved`` reads interleaved blocks.
Returns each alignment's taxa, length and sequences (at most ``limit`` blocks).
```

**输出：** 返回可 JSON 序列化的字典；文件输出返回路径。可选依赖缺失时抛出带安装提示的 BackendUnavailableError。

## `dnakit.phylogeny.dendropy_p38_functions.calculate_hky85_rate_matrix`

**作用：** HKY85 instantaneous nucleotide rate matrix Q and the branch-length correction factor, using

**处理方法：** 由模块中的显式接口校验输入并调用对应后端。

**API：**

```text
def calculate_hky85_rate_matrix(*, base_frequencies: Sequence[float]=(0.25, 0.25, 0.25, 0.25), kappa: float=1.0, rate: float=1.0, branch_length: float | None=None, verbose: bool=True) -> dict[str, Any]:
```

**输入：** 签名列出必选参数、可选参数与原生默认值；参数含义和输出如下。

```text
HKY85 instantaneous nucleotide rate matrix Q and the branch-length correction factor, using
``dendropy.model.discrete.Hky85.qmatrix`` and ``corrected_substitution_rate``.

``base_frequencies`` are ordered A, C, G, T; ``kappa`` is the transition/transversion rate ratio; ``rate`` is
the relative rate. ``corrected_substitution_rate`` is the factor that converts a branch length to expected
substitutions per site (``rate * 0.5 / f``). With ``branch_length`` the matrix of transition probabilities
over that length (``pij``) is returned too.
```

**输出：** 返回可 JSON 序列化的字典；文件输出返回路径。可选依赖缺失时抛出带安装提示的 BackendUnavailableError。

## `dnakit.phylogeny.dendropy_p38_functions.sample_stationary_nucleotide_sequence`

**作用：** Sample a root nucleotide sequence from the stationary base distribution, using

**处理方法：** 由模块中的显式接口校验输入并调用对应后端。

**API：**

```text
def sample_stationary_nucleotide_sequence(length: int, *, base_frequencies: Sequence[float]=(0.25, 0.25, 0.25, 0.25), seed: int | None=None, verbose: bool=True) -> dict[str, Any]:
```

**输入：** 签名列出必选参数、可选参数与原生默认值；参数含义和输出如下。

```text
Sample a root nucleotide sequence from the stationary base distribution, using
``dendropy.model.discrete.NucleotideCharacterEvolutionModel.stationary_sample``.

``base_frequencies`` are ordered A, C, G, T; ``seed`` makes the draw reproducible.
```

**输出：** 返回可 JSON 序列化的字典；文件输出返回路径。可选依赖缺失时抛出带安装提示的 BackendUnavailableError。

## `dnakit.phylogeny.dendropy_p38_functions.simulate_fragmented_population_sequences`

**作用：** Simulate DNA sequences from fragmented (structured) populations, using

**处理方法：** 由模块中的显式接口校验输入并调用对应后端。

**API：**

```text
def simulate_fragmented_population_sequences(species_name: str, *, div_time_gens: int, num_desc_pops: int=2, mutrate_per_site_per_generation: float=1e-07, desc_pop_size: int=10000, samples_per_pop: int=10, seq_len: int=2000, seed: int | None=None, output_path: str | None=None, output_format: str='fasta', verbose: bool=True) -> dict[str, Any]:
```

**输入：** 签名列出必选参数、可选参数与原生默认值；参数含义和输出如下。

```text
Simulate DNA sequences from fragmented (structured) populations, using
``dendropy.simulate.popgensim.FragmentedPopulations.generate_sequences``.

A population tree with ``num_desc_pops`` descendant populations of ``desc_pop_size`` diverging
``div_time_gens`` generations ago is built, a gene tree is drawn with ``samples_per_pop`` genes per
population, and ``seq_len`` sites are evolved on the mutation tree (dendropy's built-in HKY85, equal base
frequencies). The sequences are returned and, with ``output_path``, written in ``output_format`` (fasta,
nexus, nexml or phylip). ``seed`` makes the run reproducible.
```

**输出：** 返回可 JSON 序列化的字典；文件输出返回路径。可选依赖缺失时抛出带安装提示的 BackendUnavailableError。

## `dnakit.phylogeny.dendropy_p38_functions.calculate_folded_site_frequency_spectrum`

**作用：** Folded (minor-allele) site frequency spectrum of a DNA alignment, using

**处理方法：** 由模块中的显式接口校验输入并调用对应后端。

**API：**

```text
def calculate_folded_site_frequency_spectrum(sequences: Mapping[str, str], *, pad_to_unfolded_length: bool=False, verbose: bool=True) -> dict[str, Any]:
```

**输入：** 签名列出必选参数、可选参数与原生默认值；参数含义和输出如下。

```text
Folded (minor-allele) site frequency spectrum of a DNA alignment, using
``DiscreteCharacterMatrix.folded_site_frequency_spectrum``.

``sequences`` maps IDs to aligned DNA strings. Entry i counts the sites whose minor allele occurs in i
sequences; the vector has ceil(N/2) entries, or N+1 with ``pad_to_unfolded_length``, whose first entry then
counts monomorphic sites (as dendropy does).
```

**输出：** 返回可 JSON 序列化的字典；文件输出返回路径。可选依赖缺失时抛出带安装提示的 BackendUnavailableError。

## `dnakit.phylogeny.dendropy_p38_functions.calculate_population_pair_statistics`

**作用：** Two-population pairwise summary statistics, using

**处理方法：** 由模块中的显式接口校验输入并调用对应后端。

**API：**

```text
def calculate_population_pair_statistics(population1: Mapping[str, str], population2: Mapping[str, str], *, ignore_uncertain: bool=True, verbose: bool=True) -> dict[str, Any]:
```

**输入：** 签名列出必选参数、可选参数与原生默认值；参数含义和输出如下。

```text
Two-population pairwise summary statistics, using
``dendropy.calculate.popgenstat.PopulationPairSummaryStatistics`` (Hickerson et al. 2006, Wakeley 1996).

Each population maps IDs to aligned DNA strings of the same length. Returns average pairwise differences
(overall, between, within, net), segregating sites, Watterson's theta, Wakeley's psi and Tajima's D. Needs
at least two sequences per population.
```

**输出：** 返回可 JSON 序列化的字典；文件输出返回路径。可选依赖缺失时抛出带安装提示的 BackendUnavailableError。

## `dnakit.phylogeny.dendropy_p38_functions.parse_genbank_xml_records`

**作用：** Parse GenBank INSDSeq XML records offline, using dendropy's GenBank interop classes.

**处理方法：** 由模块中的显式接口校验输入并调用对应后端。

**API：**

```text
def parse_genbank_xml_records(xml_path: str, *, molecule: str='dna', include_features: bool=True, fasta_labels: str='defline', limit: int=100, verbose: bool=True) -> dict[str, Any]:
```

**输入：** 签名列出必选参数、可选参数与原生默认值；参数含义和输出如下。

```text
Parse GenBank INSDSeq XML records offline, using dendropy's GenBank interop classes.

``xml_path`` is an NCBI efetch ``rettype=gbc, retmode=xml`` file (INSDSet/INSDSeq XML, the format dendropy's
GenBank interop reads; plain ``gb`` GBSeq XML is not read by dendropy). For each record the function returns
the metadata (accession, version, organism, definition, taxonomy, dates, other sequence IDs), the sequence,
its FASTA text (``GenBankAccessionRecord.as_fasta``; ``fasta_labels`` is ``defline`` for the NCBI-style
``gi|..|gb|..|`` header or ``accession_organism`` for a composed label) and, with ``include_features``, the
feature table (key, location, qualifiers). ``molecule`` is ``dna`` or ``rna``.
```

**输出：** 返回可 JSON 序列化的字典；文件输出返回路径。可选依赖缺失时抛出带安装提示的 BackendUnavailableError。

## `dnakit.phylogeny.dendropy_p38_functions.genbank_xml_to_alignment`

**作用：** Build a nucleotide character matrix from GenBank INSDSeq XML records and write it, using

**处理方法：** 由模块中的显式接口校验输入并调用对应后端。

**API：**

```text
def genbank_xml_to_alignment(xml_path: str, output_path: str, *, output_format: str='nexml', molecule: str='dna', label_components: Sequence[str]=('accession', 'organism'), label_separator: str='_', verbose: bool=True) -> dict[str, Any]:
```

**输入：** 签名列出必选参数、可选参数与原生默认值；参数含义和输出如下。

```text
Build a nucleotide character matrix from GenBank INSDSeq XML records and write it, using
``GenBankResourceStore.generate_char_matrix``.

Taxon labels are composed from ``label_components`` (record attributes such as ``accession``, ``organism``,
``gi``) joined by ``label_separator``. The records need not be aligned: write ``fasta`` for unaligned sequences.
(dendropy 5.1's GenBank annotation helpers ``as_annotation`` and ``as_reference_annotation`` raise ``TypeError``
and are not used; metadata and feature tables are available from ``parse_genbank_xml_records``.)
```

**输出：** 返回可 JSON 序列化的字典；文件输出返回路径。可选依赖缺失时抛出带安装提示的 BackendUnavailableError。

## `dnakit.phylogeny.dendropy_p38_functions.concatenate_alignments`

**作用：** Concatenate alignments into a partitioned supermatrix, using

**处理方法：** 由模块中的显式接口校验输入并调用对应后端。

**API：**

```text
def concatenate_alignments(alignment_paths: Sequence[str], output_path: str, *, input_format: str='fasta', output_format: str='nexus', data_type: str='dna', verbose: bool=True) -> dict[str, Any]:
```

**输入：** 签名列出必选参数、可选参数与原生默认值；参数含义和输出如下。

```text
Concatenate alignments into a partitioned supermatrix, using
``dendropy.datamodel.charmatrixmodel.CharacterMatrix.concatenate`` (and ``concatenate_from_paths``).

Every alignment must contain the same taxa (labels are matched by name). The parts are recorded as
character subsets named after the input files; the ``nexus`` output keeps them as ``charset`` blocks. The
result reports each part's column range (0-based, half-open).
```

**输出：** 返回可 JSON 序列化的字典；文件输出返回路径。可选依赖缺失时抛出带安装提示的 BackendUnavailableError。

## `dnakit.phylogeny.dendropy_p38_functions.extract_alignment_columns`

**作用：** Extract alignment columns by index or range, using ``CharacterMatrix.export_character_indices`` (and

**处理方法：** 由模块中的显式接口校验输入并调用对应后端。

**API：**

```text
def extract_alignment_columns(alignment_path: str, output_path: str, *, columns: Sequence[int] | None=None, ranges: Sequence[Sequence[int]] | None=None, input_format: str='fasta', output_format: str='fasta', data_type: str='dna', verbose: bool=True) -> dict[str, Any]:
```

**输入：** 签名列出必选参数、可选参数与原生默认值；参数含义和输出如下。

```text
Extract alignment columns by index or range, using ``CharacterMatrix.export_character_indices`` (and
``export_character_subset`` for named subsets).

``columns`` lists 0-based column indices; ``ranges`` lists ``[start, end)`` intervals (0-based, half-open);
both can be given and are merged in alignment order. When the alignment is NEXUS with ``charset`` blocks,
``input_format="nexus"`` keeps them and ``ranges`` can be replaced by naming nothing: the function
extracts only the requested columns.
```

**输出：** 返回可 JSON 序列化的字典；文件输出返回路径。可选依赖缺失时抛出带安装提示的 BackendUnavailableError。

## 可执行示例

**示例代码：**

```python
from dnakit.molecular_lab.dnachisel_p38_functions import evaluate_dnachisel_design_rule
from dnakit.phylogeny.dendropy_p38_functions import calculate_hky85_rate_matrix
r = evaluate_dnachisel_design_rule('ACGT', 'SequenceLengthBounds', parameters={'min_length': 2, 'max_length': 6}, verbose=False)
print(r['evaluation']['score'])
print(calculate_hky85_rate_matrix(verbose=False)['qmatrix'][0])
```

**示例输出：**

```text
0.0
[-1.0, 0.3333333333333333, 0.3333333333333333, 0.3333333333333333]
```
