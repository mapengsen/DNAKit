# minigraph: sequence-to-graph mapping, graph construction, bubble calls and SV scripts

This group of functions in `dnakit.alignment` is a set of command-line facades over [minigraph](https://github.com/lh3/minigraph) and its companion scripts. `map_sequences_to_graph`, `graph_segment_coverage`, `build_sequence_graph` and `call_graph_bubbles` run the `minigraph` program (mapping, `--cov`, `-cxggs`, `--call`). The other functions run one sub-command of `mgutils.js` (rename FASTA, paf2bl, anno, anno2tbl, stableGaf, extractseg, merge, segfreq, genecopy), `mgutils-es6.js` (merge2vcf, addsample) or `gafcall.js` (extract, merge, eval, view, join) through the `k8` JavaScript shell. `mgutils.js joinfa` is not provided: it fails with every released k8 (it calls `Bytes.set` with a `Bytes` argument, which k8 1.2 rejects).

minigraph 0.21 and k8 1.2 are the tested versions (MIT licence; started as separate programs). Install them with `conda install -c bioconda minigraph k8` or build them from the project pages; put `minigraph`, `k8` and the scripts (`mgutils.js`, `mgutils-es6.js`, `gafcall.js`, in the `misc/` folder of the source tree) on PATH or pass `minigraph_path`, `k8_path`, `mgutils_path`, `mgutils_es6_path` or `gafcall_path`. Every function writes into a new `output_directory`; programs that print to standard output are streamed into a file there, standard error goes to `tool.log`. Options left at `None` are not passed, so the tool's defaults apply. The examples use `verbose=False`.

## Example data and examples

**Example code:**

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

**Example output:**

```text
0.21 3 3
1 ['ref', '2999', '2999']
1 1
```

**Example code:**

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

**Example output:**

```text
1
2 Other-repeat
```


## 1) `dnakit.alignment.minigraph_p11_cli_functions.add_sample_names_to_calls()`

- **Purpose:** Put sample names into the header of a merged bubble-call table (``mgutils-es6.js addsample``).
- **Method:** Replaces the ``#CHROM`` header line of ``merged_calls`` with one that lists ``samples``.
- **API:** `dnakit.alignment.minigraph_p11_cli_functions.add_sample_names_to_calls(merged_calls: 'str | Path', samples: 'Sequence[str]', output_directory: 'str | Path', *, output_name: 'str' = 'merged_named.txt', mgutils_es6_path: 'str | None' = None, k8_path: 'str | None' = None, timeout_seconds: 'float' = 3600.0, verbose: 'bool' = True) -> 'dict[str, Any]'`.
- **Input:** required: `merged_calls`, `samples`, `output_directory`; every other parameter is keyword-only and optional (see the API line; `None` or `False` means the tool's own default).
- **Output:** a dictionary with `status, backend, operation, version, command, output_directory, outputs, log_tail, provenance`; Result key: ``sample_count``.

## 2) `dnakit.alignment.minigraph_p11_cli_functions.annotate_inserted_sequences()`

- **Purpose:** Annotate short inserted sequences (repeats, low-complexity regions, gaps, segmental duplications; ``mgutils.js anno``).
- **Method:** ``bed`` lists the inserted sequences as ``contig start end length``. Evidence files (all optional): the RepeatMasker ``.out`` file (``-r``), the stretches of N reported by ``seqtk cutN``/``seqtk gap`` (``-g``), minimap2/sdust output for low-complexity regions (``-d``), eTRF output (``-e``), a PAF alignment against the reference (``-p``), dna-brnn centromere calls (``-c``), the ``gfatools bubble`` file (``-b``) and the segmental-duplication table (``-s``). ``min_feature_length`` is ``-l`` (default 30), ``min_segdup_fraction`` is ``-S`` (default 0.2) and ``annotate_mobile`` (``-m``) reports AluY and L1HS separately.
- **API:** `dnakit.alignment.minigraph_p11_cli_functions.annotate_inserted_sequences(bed: 'str | Path', output_directory: 'str | Path', *, output_name: 'str' = 'anno.tsv', repeatmasker_out: 'str | Path | None' = None, gaps: 'str | Path | None' = None, low_complexity: 'str | Path | None' = None, tandem_repeats: 'str | Path | None' = None, paf: 'str | Path | None' = None, centromere: 'str | Path | None' = None, bubbles: 'str | Path | None' = None, segdup: 'str | Path | None' = None, min_feature_length: 'int | None' = None, min_segdup_fraction: 'float | None' = None, annotate_mobile: 'bool' = False, mgutils_path: 'str | None' = None, k8_path: 'str | None' = None, timeout_seconds: 'float' = 3600.0, verbose: 'bool' = True) -> 'dict[str, Any]'`.
- **Input:** required: `bed`, `output_directory`; every other parameter is keyword-only and optional (see the API line; `None` or `False` means the tool's own default).
- **Output:** a dictionary with `status, backend, operation, version, command, output_directory, outputs, log_tail, provenance`; Result key: ``annotated`` (rows written).

## 3) `dnakit.alignment.minigraph_p11_cli_functions.blacklist_insert_alignments()`

- **Purpose:** Blacklist reference regions that inserted sequences align to (``mgutils.js paf2bl``).
- **Method:** ``paf`` holds alignments of the inserted sequences to the reference (with the ``tp`` and ``de`` tags written by minimap2). Alignments with mapping quality below 5 or length below 500 are ignored. By default regions with a divergence between ``min_divergence`` (``-d``, default 0.01) and 0.1 are written; ``identical_only`` (``-s``) writes regions with divergence up to 0.002 instead. The output is BED (reference, start, end).
- **API:** `dnakit.alignment.minigraph_p11_cli_functions.blacklist_insert_alignments(paf: 'str | Path', output_directory: 'str | Path', *, output_name: 'str' = 'blacklist.bed', min_divergence: 'float' = 0.01, identical_only: 'bool' = False, mgutils_path: 'str | None' = None, k8_path: 'str | None' = None, timeout_seconds: 'float' = 3600.0, verbose: 'bool' = True) -> 'dict[str, Any]'`.
- **Input:** required: `paf`, `output_directory`; every other parameter is keyword-only and optional (see the API line; `None` or `False` means the tool's own default).
- **Output:** a dictionary with `status, backend, operation, version, command, output_directory, outputs, log_tail, provenance`; Result key: ``regions``.

## 4) `dnakit.alignment.minigraph_p11_cli_functions.build_sequence_graph()`

- **Purpose:** Build a sequence graph by incrementally adding assemblies to a reference (minigraph ``-cxggs``).
- **Method:** ``reference`` is a FASTA file or an existing GFA graph; the ``assemblies`` are added in the order given. ``preset`` ``ggs`` (the default) is ``-x ggs``; ``None`` runs the bare ``--ggen`` mode without the preset parameters. ``min_mapq`` (``-q``), ``min_aln_length`` (``-l``), ``min_depth_aln_length`` (``-d``) and ``min_variant_length`` (``-L``) keep minigraph's defaults (5, 100000, 20000 and 50) unless set. The output is GFA.
- **API:** `dnakit.alignment.minigraph_p11_cli_functions.build_sequence_graph(reference: 'str | Path', assemblies: 'Sequence[str | Path]', output_directory: 'str | Path', *, output_name: 'str' = 'graph.gfa', preset: 'str | None' = 'ggs', base_alignment: 'bool' = True, min_mapq: 'int | None' = None, min_aln_length: 'int | None' = None, min_depth_aln_length: 'int | None' = None, min_variant_length: 'int | None' = None, kmer_size: 'int | None' = None, window_size: 'int | None' = None, threads: 'int' = 4, minigraph_path: 'str | None' = None, timeout_seconds: 'float' = 3600.0, verbose: 'bool' = True) -> 'dict[str, Any]'`.
- **Input:** required: `reference`, `assemblies`, `output_directory`; every other parameter is keyword-only and optional (see the API line; `None` or `False` means the tool's own default).
- **Output:** a dictionary with `status, backend, operation, version, command, output_directory, outputs, log_tail, provenance`; Result keys: ``segments``, ``links``, ``paths``.

## 5) `dnakit.alignment.minigraph_p11_cli_functions.call_graph_bubbles()`

- **Purpose:** Call the path an assembly takes through every bubble of a graph, one BED row per bubble (minigraph ``--call``).
- **Method:** This is the per-sample step of the minigraph pangenome pipeline (``minigraph -cxasm --call graph.gfa sample.fa``). The tool's own defaults for the alignment length thresholds (``-l`` 100000, ``-d`` 20000) are kept unless set; the minigraph documentation uses ``-l10000 -d10000`` for assemblies with short contigs.
- **API:** `dnakit.alignment.minigraph_p11_cli_functions.call_graph_bubbles(graph: 'str | Path', assembly: 'str | Path', output_directory: 'str | Path', *, output_name: 'str' = 'calls.bed', preset: 'str | None' = 'asm', base_alignment: 'bool' = True, min_mapq: 'int | None' = None, min_aln_length: 'int | None' = None, min_depth_aln_length: 'int | None' = None, min_variant_length: 'int | None' = None, kmer_size: 'int | None' = None, window_size: 'int | None' = None, threads: 'int' = 4, minigraph_path: 'str | None' = None, timeout_seconds: 'float' = 3600.0, verbose: 'bool' = True) -> 'dict[str, Any]'`.
- **Input:** required: `graph`, `assembly`, `output_directory`; every other parameter is keyword-only and optional (see the API line; `None` or `False` means the tool's own default).
- **Output:** a dictionary with `status, backend, operation, version, command, output_directory, outputs, log_tail, provenance`; Result keys: ``bubbles`` (rows), ``rows`` (first 1000 rows split into columns).

## 6) `dnakit.alignment.minigraph_p11_cli_functions.convert_calls_to_vcf()`

- **Purpose:** Convert the merged bubble calls to VCF (``mgutils-es6.js merge2vcf``).
- **Method:** ``reference_index`` is ``-r`` (which sample column is the reference, default 0), ``max_alleles`` is ``-a`` (default 15) and ``samples`` the sample names (``-s``; one per call column).
- **API:** `dnakit.alignment.minigraph_p11_cli_functions.convert_calls_to_vcf(merged_calls: 'str | Path', output_directory: 'str | Path', *, output_name: 'str' = 'calls.vcf', reference_index: 'int | None' = None, max_alleles: 'int | None' = None, samples: 'Sequence[str] | None' = None, mgutils_es6_path: 'str | None' = None, k8_path: 'str | None' = None, timeout_seconds: 'float' = 3600.0, verbose: 'bool' = True) -> 'dict[str, Any]'`.
- **Input:** required: `merged_calls`, `output_directory`; every other parameter is keyword-only and optional (see the API line; `None` or `False` means the tool's own default).
- **Output:** a dictionary with `status, backend, operation, version, command, output_directory, outputs, log_tail, provenance`; Result key: ``variants`` (VCF records).

## 7) `dnakit.alignment.minigraph_p11_cli_functions.convert_to_stable_gaf()`

- **Purpose:** Convert graph alignments (GAF) to stable coordinates on the linear sequences of the graph (``mgutils.js stableGaf``).
- **Method:** ``graph`` is the minigraph GFA with the ``SN``, ``SO`` and ``SR`` segment tags; ``alignments`` is the GAF written by minigraph. Alignments inside one reference-rank-0 segment path are rewritten in the coordinates of that sequence.
- **API:** `dnakit.alignment.minigraph_p11_cli_functions.convert_to_stable_gaf(graph: 'str | Path', alignments: 'str | Path', output_directory: 'str | Path', *, output_name: 'str' = 'stable.gaf', mgutils_path: 'str | None' = None, k8_path: 'str | None' = None, timeout_seconds: 'float' = 3600.0, verbose: 'bool' = True) -> 'dict[str, Any]'`.
- **Input:** required: `graph`, `alignments`, `output_directory`; every other parameter is keyword-only and optional (see the API line; `None` or `False` means the tool's own default).
- **Output:** a dictionary with `status, backend, operation, version, command, output_directory, outputs, log_tail, provenance`; Result key: ``records``.

## 8) `dnakit.alignment.minigraph_p11_cli_functions.count_gene_copies()`

- **Purpose:** Count the copies of each gene in an assembly from graph alignments of the gene sequences (``mgutils.js genecopy``).
- **Method:** ``alignments`` is the GAF (with ``cg`` tags) of gene sequences mapped to the assembly graph, ``gene_bed`` the BED (name in column 4, strand in column 6) of the genes on the source reference. ``min_coverage`` is ``-c`` (default 0.8) and ``min_relative_coverage`` is ``-r`` (default 0.85).
- **API:** `dnakit.alignment.minigraph_p11_cli_functions.count_gene_copies(alignments: 'str | Path', gene_bed: 'str | Path', output_directory: 'str | Path', *, output_name: 'str' = 'gene_copy.tsv', min_coverage: 'float | None' = None, min_relative_coverage: 'float | None' = None, mgutils_path: 'str | None' = None, k8_path: 'str | None' = None, timeout_seconds: 'float' = 3600.0, verbose: 'bool' = True) -> 'dict[str, Any]'`.
- **Input:** required: `alignments`, `gene_bed`, `output_directory`; every other parameter is keyword-only and optional (see the API line; `None` or `False` means the tool's own default).
- **Output:** a dictionary with `status, backend, operation, version, command, output_directory, outputs, log_tail, provenance`; Result key: ``rows`` (output lines, first 1000).

## 9) `dnakit.alignment.minigraph_p11_cli_functions.evaluate_sv_calls()`

- **Purpose:** Evaluate structural-variant calls against a truth set, or compare several call sets (``gafcall.js eval``).
- **Method:** With two VCF files, the first is the truth ("base") set and the second the call ("test") set; the result has the recall (``recall``: truth SVs missed and found) and precision-side counts (``precision``). With more than two files the result is the pairwise concordance matrix (``concordance``). Optional parameters (gafcall defaults in brackets): ``confident_regions`` BED (``-b``), ``min_sv_length`` (``-l``, 100), ``window_size`` (``-w``, 500), ``read_length_ratio`` (``-r``, 0.8), ``min_length_ratio`` (``-m``, 0.6), ``min_allele_frequency`` (``-v``, 0), ``ignore_filter`` (``-F``), ``check_genotype`` (``-G``) and ``print_errors`` (``-e``). The raw report is in ``outputs['output']``.
- **API:** `dnakit.alignment.minigraph_p11_cli_functions.evaluate_sv_calls(vcf_files: 'Sequence[str | Path]', output_directory: 'str | Path', *, output_name: 'str' = 'eval.txt', confident_regions: 'str | Path | None' = None, min_sv_length: 'int | None' = None, window_size: 'int | None' = None, read_length_ratio: 'float | None' = None, min_length_ratio: 'float | None' = None, min_allele_frequency: 'float | None' = None, ignore_filter: 'bool' = False, check_genotype: 'bool' = False, print_errors: 'bool' = False, gafcall_path: 'str | None' = None, k8_path: 'str | None' = None, timeout_seconds: 'float' = 3600.0, verbose: 'bool' = True) -> 'dict[str, Any]'`.
- **Input:** required: `vcf_files`, `output_directory`; every other parameter is keyword-only and optional (see the API line; `None` or `False` means the tool's own default).
- **Output:** a dictionary with `status, backend, operation, version, command, output_directory, outputs, log_tail, provenance`; Result keys: see Method.

## 10) `dnakit.alignment.minigraph_p11_cli_functions.extract_gaf_segment_path()`

- **Purpose:** Extract the contig interval between two graph segments from GAF alignments (``mgutils.js extractseg``).
- **Method:** ``first_segment`` and ``second_segment`` are segment names (for example ``s12``). For every alignment file the interval of each contig that lies between the two segments is written as BED (contig, start, end, ``*``, 0, strand). Alignments shorter than ``min_alignment_length`` (``-l``, default 100000) are ignored; ``from_end`` (``-e``) takes the interval from the end of the first segment to the start of the second.
- **API:** `dnakit.alignment.minigraph_p11_cli_functions.extract_gaf_segment_path(first_segment: 'str', second_segment: 'str', alignments: 'Sequence[str | Path]', output_directory: 'str | Path', *, output_name: 'str' = 'segment_path.bed', min_alignment_length: 'int' = 100000, from_end: 'bool' = False, mgutils_path: 'str | None' = None, k8_path: 'str | None' = None, timeout_seconds: 'float' = 3600.0, verbose: 'bool' = True) -> 'dict[str, Any]'`.
- **Input:** required: `first_segment`, `second_segment`, `alignments`, `output_directory`; every other parameter is keyword-only and optional (see the API line; `None` or `False` means the tool's own default).
- **Output:** a dictionary with `status, backend, operation, version, command, output_directory, outputs, log_tail, provenance`; Result key: ``intervals``.

## 11) `dnakit.alignment.minigraph_p11_cli_functions.extract_long_indels()`

- **Purpose:** Extract long INDELs and breakpoints from stable-coordinate GAF alignments of reads (``gafcall.js extract``).
- **Method:** ``stable_gaf`` comes from ``convert_to_stable_gaf``. Optional parameters (gafcall defaults in brackets): ``sample_name`` (``-n``, foo), ``min_indel_length`` (``-l``, 100), ``max_indels_per_10kb`` (``-c``, 3), ``min_mapq`` (``-q``, 5), ``min_mapq_ends`` (``-Q``, 30), ``min_alignment_length_ends`` (``-e``, 2000), ``min_alignment_length_middle`` (``-m``, 50), ``polya_penalty`` (``-a``, 5) and a centromere BED (``-b``).
- **API:** `dnakit.alignment.minigraph_p11_cli_functions.extract_long_indels(stable_gaf: 'str | Path', output_directory: 'str | Path', *, output_name: 'str' = 'extract.tsv', sample_name: 'str | None' = None, min_indel_length: 'int | None' = None, max_indels_per_10kb: 'int | None' = None, min_mapq: 'int | None' = None, min_mapq_ends: 'int | None' = None, min_alignment_length_ends: 'int | None' = None, min_alignment_length_middle: 'int | None' = None, polya_penalty: 'int | None' = None, centromere_bed: 'str | Path | None' = None, gafcall_path: 'str | None' = None, k8_path: 'str | None' = None, timeout_seconds: 'float' = 3600.0, verbose: 'bool' = True) -> 'dict[str, Any]'`.
- **Input:** required: `stable_gaf`, `output_directory`; every other parameter is keyword-only and optional (see the API line; `None` or `False` means the tool's own default).
- **Output:** a dictionary with `status, backend, operation, version, command, output_directory, outputs, log_tail, provenance`; Result key: ``records`` (evidence lines).

## 12) `dnakit.alignment.minigraph_p11_cli_functions.filter_sv_extract_by_type()`

- **Purpose:** Keep the records of one SV-extract file whose read also carries the same SV type in another file (``gafcall.js join``).
- **Method:** A record of ``extract`` is kept when its read name occurs in ``filter_extract`` with a matching type (insertion or duplication, deletion, ...) or when the type is undefined.
- **API:** `dnakit.alignment.minigraph_p11_cli_functions.filter_sv_extract_by_type(filter_extract: 'str | Path', extract: 'str | Path', output_directory: 'str | Path', *, output_name: 'str' = 'joined.tsv', gafcall_path: 'str | None' = None, k8_path: 'str | None' = None, timeout_seconds: 'float' = 3600.0, verbose: 'bool' = True) -> 'dict[str, Any]'`.
- **Input:** required: `filter_extract`, `extract`, `output_directory`; every other parameter is keyword-only and optional (see the API line; `None` or `False` means the tool's own default).
- **Output:** a dictionary with `status, backend, operation, version, command, output_directory, outputs, log_tail, provenance`; Result key: ``records`` (lines kept).

## 13) `dnakit.alignment.minigraph_p11_cli_functions.graph_segment_coverage()`

- **Purpose:** Remap reads (or assemblies) to a graph and report the coverage of every graph segment (minigraph ``--cov``).
- **Method:** The output is the input graph in GFA format with a ``dc:f:`` (mean depth) tag on every covered segment. ``min_aligned_length`` and ``min_mapq`` are ``--min-cov-blen`` and ``--min-cov-mapq`` (alignments below them do not count).
- **API:** `dnakit.alignment.minigraph_p11_cli_functions.graph_segment_coverage(graph: 'str | Path', reads: 'Sequence[str | Path]', output_directory: 'str | Path', *, output_name: 'str' = 'coverage.gfa', preset: 'str | None' = None, base_alignment: 'bool' = True, min_aligned_length: 'int | None' = None, min_mapq: 'int | None' = None, kmer_size: 'int | None' = None, window_size: 'int | None' = None, threads: 'int' = 4, minigraph_path: 'str | None' = None, timeout_seconds: 'float' = 3600.0, verbose: 'bool' = True) -> 'dict[str, Any]'`.
- **Input:** required: `graph`, `reads`, `output_directory`; every other parameter is keyword-only and optional (see the API line; `None` or `False` means the tool's own default).
- **Output:** a dictionary with `status, backend, operation, version, command, output_directory, outputs, log_tail, provenance`; Result keys: ``segments``, ``covered_segments`` (segments carrying the depth tag), ``mean_depths`` (segment name to depth, at most 1000 segments).

## 14) `dnakit.alignment.minigraph_p11_cli_functions.map_sequences_to_graph()`

- **Purpose:** Map DNA sequences (reads or assemblies) to a sequence graph, or to a linear reference (minigraph).
- **Method:** ``graph`` is a GFA file (or a FASTA file for a linear target); ``queries`` are FASTA/FASTQ files, optionally gzipped. ``preset`` is ``lr`` (noisy long reads, the default), ``asm`` (assembly-to-graph), ``sr`` (short reads) or ``ggs``. ``base_alignment`` is minigraph's ``-c`` (base-level alignment; recommended by the authors). The output is GAF (PAF when ``graph`` is linear) written to ``<output_directory>/<output_name>``; with ``vertex_coordinates`` the path coordinates are vertex coordinates.
- **API:** `dnakit.alignment.minigraph_p11_cli_functions.map_sequences_to_graph(graph: 'str | Path', queries: 'Sequence[str | Path]', output_directory: 'str | Path', *, output_name: 'str' = 'mapping.gaf', preset: 'str | None' = None, base_alignment: 'bool' = False, kmer_size: 'int | None' = None, window_size: 'int | None' = None, repetitive_fraction: 'float | None' = None, occurrence_range: 'int | Sequence[int] | None' = None, divergence: 'float | None' = None, chain_gap: 'int | None' = None, bandwidth: 'int | Sequence[int] | None' = None, min_minimizers: 'int | Sequence[int] | None' = None, min_chain_score: 'int | Sequence[int] | None' = None, min_secondary_ratio: 'float | None' = None, max_secondary: 'int | None' = None, skip_self_diagonal: 'bool' = False, linear_chains: 'bool' = False, vertex_coordinates: 'bool' = False, threads: 'int' = 4, minigraph_path: 'str | None' = None, timeout_seconds: 'float' = 3600.0, verbose: 'bool' = True) -> 'dict[str, Any]'`.
- **Input:** required: `graph`, `queries`, `output_directory`; every other parameter is keyword-only and optional (see the API line; `None` or `False` means the tool's own default).
- **Output:** a dictionary with `status, backend, operation, version, command, output_directory, outputs, log_tail, provenance`; Result keys: ``records`` (number of alignment lines), ``mapped_queries`` (number of distinct query names), ``outputs``.

## 15) `dnakit.alignment.minigraph_p11_cli_functions.merge_bubble_calls()`

- **Purpose:** Merge the per-sample ``--call`` BED files into one multi-sample table (``paste *.bed | mgutils.js merge``).
- **Method:** ``call_files`` are the BED files written by ``call_graph_bubbles`` for the same graph (same bubbles in the same order), one per sample. ``samples`` are the sample names for the header (``-s``); ``annotation`` is a bubble annotation file whose twelfth column is added as ``ANNO`` (``-a``). The output has a VCF-like header, one row per bubble and one genotype column per sample.
- **API:** `dnakit.alignment.minigraph_p11_cli_functions.merge_bubble_calls(call_files: 'Sequence[str | Path]', output_directory: 'str | Path', *, output_name: 'str' = 'merged.txt', samples: 'Sequence[str] | None' = None, annotation: 'str | Path | None' = None, mgutils_path: 'str | None' = None, k8_path: 'str | None' = None, timeout_seconds: 'float' = 3600.0, verbose: 'bool' = True) -> 'dict[str, Any]'`.
- **Input:** required: `call_files`, `output_directory`; every other parameter is keyword-only and optional (see the API line; `None` or `False` means the tool's own default).
- **Output:** a dictionary with `status, backend, operation, version, command, output_directory, outputs, log_tail, provenance`; Result keys: ``bubbles`` (rows), ``sample_count``.

## 16) `dnakit.alignment.minigraph_p11_cli_functions.merge_sv_evidence()`

- **Purpose:** Merge read-level INDEL and breakpoint evidence into structural-variant calls (``gafcall.js merge``).
- **Method:** ``extract_files`` are outputs of ``extract_long_indels`` (several files, for example of different read groups of the same sample, are concatenated and sorted by contig and position first, as ``sort -k1,1 -k2,2n`` does). Optional parameters (gafcall defaults in brackets): ``min_read_count`` (``-c``, 4), ``min_read_count_per_strand`` (``-s``, 2), ``window_size`` (``-w``, 100), ``max_length_difference`` (``-d``, 0.05), ``min_centromere_distance`` (``-e``, 500000), ``min_retrotransposon_length`` (``-r``, 0 disables), ``min_retrotransposon_reads`` (``-R``, 1), ``max_alleles`` (``-A``, 100) and ``max_reads_compared`` (``-C``, 500).
- **API:** `dnakit.alignment.minigraph_p11_cli_functions.merge_sv_evidence(extract_files: 'Sequence[str | Path]', output_directory: 'str | Path', *, output_name: 'str' = 'calls.tsv', min_read_count: 'int | None' = None, min_read_count_per_strand: 'int | None' = None, window_size: 'int | None' = None, max_length_difference: 'float | None' = None, min_centromere_distance: 'int | None' = None, min_retrotransposon_length: 'int | None' = None, min_retrotransposon_reads: 'int | None' = None, max_alleles: 'int | None' = None, max_reads_compared: 'int | None' = None, gafcall_path: 'str | None' = None, k8_path: 'str | None' = None, timeout_seconds: 'float' = 3600.0, verbose: 'bool' = True) -> 'dict[str, Any]'`.
- **Input:** required: `extract_files`, `output_directory`; every other parameter is keyword-only and optional (see the API line; `None` or `False` means the tool's own default).
- **Output:** a dictionary with `status, backend, operation, version, command, output_directory, outputs, log_tail, provenance`; Result key: ``calls`` (rows).

## 17) `dnakit.alignment.minigraph_p11_cli_functions.rename_fasta_sequences()`

- **Purpose:** Rename FASTA sequences to ``<prefix><delimiter><name>`` (the PanSN convention; ``mgutils.js renamefa``).
- **Method:** An existing ``<something>#`` prefix of the sequence name is removed first. ``delimiter`` is mgutils' ``-d`` (default ``#``).
- **API:** `dnakit.alignment.minigraph_p11_cli_functions.rename_fasta_sequences(fasta: 'str | Path', prefix: 'str', output_directory: 'str | Path', *, output_name: 'str' = 'renamed.fa', delimiter: 'str' = '#', mgutils_path: 'str | None' = None, k8_path: 'str | None' = None, timeout_seconds: 'float' = 3600.0, verbose: 'bool' = True) -> 'dict[str, Any]'`.
- **Input:** required: `fasta`, `prefix`, `output_directory`; every other parameter is keyword-only and optional (see the API line; `None` or `False` means the tool's own default).
- **Output:** a dictionary with `status, backend, operation, version, command, output_directory, outputs, log_tail, provenance`; Result key: ``sequences`` (number of FASTA records written).

## 18) `dnakit.alignment.minigraph_p11_cli_functions.segment_allele_frequency()`

- **Purpose:** Compute the allele frequency of every graph segment from merged bubble calls (``mgutils.js segfreq``).
- **Method:** ``segment_bed`` is the ``gfatools gfa2bed`` table of the graph, ``merged_calls`` the output of ``merge_bubble_calls`` and ``bubble_bed`` an optional ``gfatools bubble`` file. ``min_frequency`` is ``-f`` (default 0.05).
- **API:** `dnakit.alignment.minigraph_p11_cli_functions.segment_allele_frequency(segment_bed: 'str | Path', merged_calls: 'str | Path', output_directory: 'str | Path', *, output_name: 'str' = 'segfreq.tsv', bubble_bed: 'str | Path | None' = None, min_frequency: 'float | None' = None, mgutils_path: 'str | None' = None, k8_path: 'str | None' = None, timeout_seconds: 'float' = 3600.0, verbose: 'bool' = True) -> 'dict[str, Any]'`.
- **Input:** required: `segment_bed`, `merged_calls`, `output_directory`; every other parameter is keyword-only and optional (see the API line; `None` or `False` means the tool's own default).
- **Output:** a dictionary with `status, backend, operation, version, command, output_directory, outputs, log_tail, provenance`; Result key: ``segments`` (rows).

## 19) `dnakit.alignment.minigraph_p11_cli_functions.summarize_insert_annotation()`

- **Purpose:** Summarise the output of ``annotate_inserted_sequences`` as a table of repeat classes (``mgutils.js anno2tbl``).
- **Method:** Gap rows and unplaced contigs (``chrUn``, ``_random``) are ignored. Each output row is the class key, its label and nine counts (inserts with 0, 1, 2 and 3-or-more annotated features, the summed insert length for each, and the summed annotated length).
- **API:** `dnakit.alignment.minigraph_p11_cli_functions.summarize_insert_annotation(annotation: 'str | Path', output_directory: 'str | Path', *, output_name: 'str' = 'anno_table.tsv', mgutils_path: 'str | None' = None, k8_path: 'str | None' = None, timeout_seconds: 'float' = 3600.0, verbose: 'bool' = True) -> 'dict[str, Any]'`.
- **Input:** required: `annotation`, `output_directory`; every other parameter is keyword-only and optional (see the API line; `None` or `False` means the tool's own default).
- **Output:** a dictionary with `status, backend, operation, version, command, output_directory, outputs, log_tail, provenance`; Result key: ``classes`` (rows, first 1000).

## 20) `dnakit.alignment.minigraph_p11_cli_functions.view_sv_calls()`

- **Purpose:** Re-format the structural variants of a VCF as a table with their long-INDEL counts (``gafcall.js view``).
- **Method:** One row per SV: contig, position, orientation, second contig, second position, SV type and SV length. With ``count_long`` (``-C``) the output is instead the number of SVs of at least 20 kb, 100 kb and 1 Mb and the number of translocations. ``min_length`` is ``-l`` (default 100), ``regions`` a BED of regions to keep (``-b``), and ``ignore_filter`` / ``check_genotype`` are ``-F`` / ``-G``.
- **API:** `dnakit.alignment.minigraph_p11_cli_functions.view_sv_calls(vcf: 'str | Path', output_directory: 'str | Path', *, output_name: 'str' = 'sv.tsv', min_length: 'int | None' = None, regions: 'str | Path | None' = None, ignore_filter: 'bool' = False, check_genotype: 'bool' = False, count_long: 'bool' = False, gafcall_path: 'str | None' = None, k8_path: 'str | None' = None, timeout_seconds: 'float' = 3600.0, verbose: 'bool' = True) -> 'dict[str, Any]'`.
- **Input:** required: `vcf`, `output_directory`; every other parameter is keyword-only and optional (see the API line; `None` or `False` means the tool's own default).
- **Output:** a dictionary with `status, backend, operation, version, command, output_directory, outputs, log_tail, provenance`; Result key: ``records`` (output lines).

