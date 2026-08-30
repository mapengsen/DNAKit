# API reference

This page lists the public namespaces that are implemented, importable, and covered by tests in the current source code. All coordinates are 0-based half-open intervals unless explicitly stated in the format conversion section.

If you want to view the roles, inputs, examples, outputs and status of each of the 185 functions in the current trace matrix, please enter from the [Function Module Index](features/index.md).

## Top-level and core objects {#core-objects}

Ordinary users only need to use `DNA` from the top level of `dnakit`: one or more sequences, IDs, features, metadata and topology are determined by the input and optional parameters of `DNA(...)`, and the return type remains unchanged. `DNASequence`, `DNARecord`, `DNASet` are reserved for high-level internal models and legacy code compatible interfaces.

::: dnakit

Core results `MetricResult` hold value, unit, method, algorithm version, parameters, conditions, provenance, and issues; domain results can use a more specialized immutable schema.

## I/O and metadata

::: dnakit.io

Main borders:

- Normal reading uses `read(..., mode="dna")` and returns `DNA`; large files use `read(..., mode="stream")` and returns `RecordSource` for single consumption; `read_one()` and `read_set()` are only compatible;
- `ReadConfig.max_sequence_symbols/max_input_bytes/max_json_depth/max_json_nodes` and `WriteConfig.max_output_bytes` place hard caps on logging I/O and CSV/TSV/JSON/JSONL embedded JSON;
- FASTA/FASTQ/CSV/TSV/JSON/JSONL and gzip are handled by unified `read`/`write`;
- GenBank is an explicit subset of commonly used fields and does not claim full INSDC compliance;
- GFF3, BED3–6 and AGP 2.1 use separate strict codecs;
- `FastaIndex` and `FastqIndex` only support uncompressed local files and detect stale indexes with source size, mtime and SHA-256; FASTQ limits strict four-line records;
- `export_table()`/`read_table()` reads and writes CSV/TSV/JSON/Parquet with explicit `TableSchema` bounds; reading and writing are limited to row/column/cell/input/output/decoded bytes respectively. Parquet relies on `io` extra and only accepts compression methods in the `ParquetCompression` whitelist.

CSV/TSV uses the literal `\N` to represent null by default, and an empty string is still a string; this convention can be changed by setting `null_value`/`missing_values` in pairs.

```python
from pathlib import Path
from tempfile import TemporaryDirectory

from dnakit.io import TableSchema, export_table, iter_chunks, read_table

chunks = tuple(iter_chunks(range(5), chunk_size=2))
assert chunks == ((0, 1), (2, 3), (4,))

schema = TableSchema(
    ("id", "score"),
    column_types={"id": "string", "score": "number"},
    nullable=(),
)
with TemporaryDirectory() as directory:
    path = Path(directory) / "scores.json"
    export_table([{"id": "a", "score": 0.75}], path, format="json", schema=schema)
    loaded = read_table(path, format="json", schema=schema)
    assert loaded.rows[0]["score"] == 0.75
```

## Reference genome download

::: dnakit.references

This module parses assembly accessions or species names through NCBI Datasets v2, streams genomic FASTA, and retains assembly report, MD5 checksum download provenance. For the complete download process, see [Download](features/15_download.md).

## Public database query

::: dnakit.search

`SearchConfig` Set hard caps on timeouts, response bytes, and number of records. The results are unified into immutable `QueryResult`, containing masked URLs, origin fields, pagination information, and provenance. Ensembl and UCSC regional entrances uniformly receive 0-based half-open coordinates; NCBI BLAST submits asynchronously by default and does not poll implicitly.

For complete query examples and 29 capability statuses, see [Data Query ](features/22_database_query.md).

## Public data download

::: dnakit.download

The new download is first temporarily stored in the target file system, MD5/SHA-256 is calculated and the manifest is written, and then the transaction is installed; if multiple files fail, it will be rolled back. ENA/ENCODE is verified when there is a provider MD5. `build_index()` Only executes explicit executable paths, `metadata()` can export query results to JSON/JSONL/CSV/TSV/XML. For a complete download example and 24 capability statuses, see [Download](features/15_download.md).

## Standardization and Verification {#standardization}

::: dnakit.standardize

`normalize()` retains IUPAC ambiguous bases by default, removes `U` and other non-DNA characters, and returns the original input snapshot, each step, and the original position of each modification. `normalize_gaps()` Just press `GapNormalizationConfig` to convert uppercase consecutive N-runs to a known length `Gap`; short Ns, explicit gaps, and ring origin boundaries are not guess-merged. `sequence_from_agp()` Assemble sequences from strict AGP 2.1 document/entries and caller component tables; reject missing components, discontinuous, out-of-bounds, gapped/circular components and unsupported orientations.

`validate(value, config=...)` is the only recommended legal entry: single record `DNA` returns a single report, multi-record `DNA` returns a collection report; as long as one record or a collection rule fails in the collection, `is_valid` becomes `False`. The old `validate_set()` only remains compatible.

```python
from dnakit import DNA
from dnakit.standardize import GapNormalizationConfig, normalize_gaps

result = normalize_gaps(
    DNA("ACNNNNTG", alphabet="iupac"),
    config=GapNormalizationConfig(min_run_length=4),
)
assert result.sequence.is_gapped
assert result.changes[0].symbol_interval.start == 2
```

## Sequence operations {#sequence-operations}

:::dnakit.ops

Ordinary users always call the unsuffixed names: `insert()`, `delete()`, `substitute()`, `mask()`, `trim()`, `reverse_complement()`, `rotate()` and `canonical_origin()`. Passing in `DNA` syncs the feature/base-by-base annotation and returns a new `DNA`. The original `*_record()` name is reserved only for high-level interfaces that require detailed change auditing.

```python
from dnakit import DNA
from dnakit.ops import reverse_complement, rotate

linear = DNA("AACG")
circular = DNA("AACG", topology="circular")
assert reverse_complement(linear).symbols == "CGTT"
assert rotate(circular, 2).symbols == "CGAA"
```

## OPS-010 Sequence Segmentation {#sequence-chunking}

::: dnakit.chunking

`ChunkingConfig()` uses fixed-length 1024 bp, non-overlapping, `train` tags by default. Passable
`strategy` Select `fixed`, `sliding`, `random`, `multiscale`, or `curriculum`:

```python
from dnakit import ChunkingConfig, iter_fasta_chunks

config = ChunkingConfig(strategy="sliding", length=1024, step=512)
for chunk in iter_fasta_chunks("input.fa", config=config):
    record = chunk.to_record()
    print(record.id, chunk.source_start, chunk.source_end, chunk.split)
```

After passing in `bed="regions.bed"`, BED uses 0-based half-open interval; column 4 is used as
`train`, `valid`, `test` and other split tags. When BED is not transmitted, each record in FASTA is
`train` processing. `LengthCurriculum((1024, 4096, 16384))` Only describes the training phase, calling
`to_config()` Then the corresponding multi-stage segmentation configuration can be generated.

## Descriptor {#descriptors}

::: dnakit.descriptors

Descriptors cover length, composition, GC/AT, skew, CpG, k-mer, Shannon entropy, homopolymer, window, codon, linguistic complexity and exact tandem-repeat union coverage. `all_descriptors()` Pressing `descriptor_schema_v1` returns a fixed 240 fields at a time, and logs the non-computable reason for each `None` value; DNAKit does not have DiProDB values ​​built in, and the last 60 items are only calculated if the caller explicitly loads a 15×16 table to which it is entitled. [See the dedicated page ](features/05_all_descriptors.md) for complete fields, formulas, units and sources.

## Patterns and comments {#patterns}

::: dnakit.patterns

The pattern module provides exact/IUPAC/regex/PWM, six reading frame ORFs, start/stop, fixed promoter patterns, caller-supplied TF PWM, restriction enzymes, PAM/guide, CpG island, reverse complementary palindromes, inverted/tandem repeats, STRs, and low complexity regions. It only returns sequence patterns and does not predict activity, binding strength, or editing efficiency.

## Thermodynamics {#thermodynamics}

::: dnakit.thermodynamics

Internal model reproduces versioned SantaLucia 1998 parameters for linear, no Gap, standard A/C/G/T, complete complement, 2–60 nt, and Na⁺+K⁺ total unit salt salts. `ThermodynamicConditions` will be saved with the results. `duplex_stability()` The default `backend="native"` uses this internal full complementation model; canonical mismatch/dangling heterodimer can be handled when the user provides an adapter with an explicit `ntthal_path` and selects `backend="primer3-cli"`. Modified and user-preset alignments are not supported.

`optical_properties()` Add double-chain average/explicit hypochromicity, nmol/mass per OD260 and explicit modification correction to the original single-chain ε260; `concentration_from_a260()` and `convert_oligo_quantity()` complete Beer–Lambert and unit safe concentration/substance mass/mass conversion. `binding_equilibrium()`, `theoretical_melting_curve()`, `terminal_stability()` and `cosolvent_tm_correction()` provide ideal two-state equilibria, theoretical curves, terminal windows and explicit empirical corrections.

`Primer3CLIAdapter` Only execute `oligotm`/`ntthal` paths explicitly given by the user. It does not import Python bindings, search from `PATH`, or automatically install or download; the unified API checks the returned results for sequence, complete conditions, max-loop, and structure options:

```python
from dnakit import DNASequence
from dnakit.thermodynamics import Primer3CLIAdapter, duplex_stability

adapter = Primer3CLIAdapter(
    ntthal_path="/opt/primer3/src/ntthal",
    thermodynamic_parameters_path="/opt/primer3/src/primer3_config",
)
result = adapter.self_dimer(DNASequence("GCGTTTTTCGC"))
print(result.delta_g_kcal_per_mol)
mismatch = duplex_stability(
    DNASequence("GTGCAT"),
    DNASequence("ATGCAA"),
    backend="primer3-cli",
    adapter=adapter,
)
assert not mismatch.fully_complementary
assert mismatch.provenance.implementation.label == "adapter"
```

## Secondary structure

::: dnakit.secondary_structure

`analyze_dot_bracket()` and probabilistically derived indicators do not rely on external backends. `probe_nupack()` is only passively positioned; `NupackAdapter` is only explicitly executed by the caller after the user has obtained the applicable license and installed it independently. The project is not automatically installed, downloaded or called silently. The current project environment does not have NUPACK, so there are no true NUPACK numerical differential conclusions.

## Three-dimensional structure

::: dnakit.structure3d

PDB native coordinate analysis only reports geometry that is computable from explicit coordinates; standard local 12-argument 3DNA output is read via `read_3dna_bp_step()`. When necessary models are missing, grooves, charges, mechanical modulus, etc. are retained as conditional functions and do not guess values ​​from ordinary sequences.

## Sequence characterization, fingerprinting and preprocessing {#fingerprints}

::: dnakit.fingerprints

Provides integer/one-hot encoding, normal or Canonical exact k-mer features, MinHash/FracMinHash k-mer Sketch, as well as motif, restriction enzyme, GC space, repeat, coding, internal thermodynamics, hybrid and multi-scale fingerprints. `FeaturePreprocessor` Fits only training data, supports missing values, standard/min-max/L1/L2 and low variance filtering.

The thermodynamic fingerprint uses a fixed 16-dimensional v2 schema: internal Tm/ΔH/ΔS/ΔG plus available/found/Tm/ΔG of hairpin/self-dimer/heterodimer. Structural items are only executed when the caller explicitly passes in the adapter; otherwise, they are processed according to the `zero`, `sentinel` or `error` missing policy, and the available bit prevents the fill value from being misread as a real calculation.

## DNA basic model characterization {#model-representations}

::: dnakit.representations

`extract_representations()` Returns a read-only float32 for each `DNARecord` by selected model
rep. LucaOne is used by default; downloaded to the running directory when checkpoint is missing
`ckpt/lucaone-gene-step36-8m/`, reused when there is a complete checkpoint. LucaOne needs to be loaded
checkpoint comes with code, so standard backends still require `allow_remote_code=True` to be set explicitly.
See all 11 models, checkpoints, dependencies and remote code boundaries
[Neural Network Representation](features/08_fingerprints.md#neural-representations); rep → k-means see
[Neural network clustering](features/10_clustering.md#data-027-neural-clustering).

## Similarity, search and comparison {#similarity-alignment}

::: dnakit.similarity

::: dnakit.alignment

`approximate_search()` supports bounded mismatch/indel; `align_pairwise()` supports global, local and double-ended free-end `semi_global`, optional linear or three-state affine gap, and returns identity/coverage. IUPAC compares as literal symbol; DP is subject to `max_cells`. MinHash indexes are deterministic exact scans of memory sketches and are not equivalent to ANN databases.

`DashingAdapter` is the opt-in external scientific computing adapter. It only accepts executable paths explicitly provided by the caller, not automatically selected from `PATH` or project third-party directories. `matrix()` Can perform Dashing exact k-mer set or HLL-sketch Jaccard; `top_k()` Stable sorting by score and original index on the same verified matrix. The memory sequence must be linear, Gap-free and of length no less than k. An explicit FASTA/FASTQ path can also be passed in. Fixed `dist` command and flag whitelisting without shell and limiting number of items, input/output/capture bytes, sketch memory, threads, and timeouts; results in logging input/raw output SHA-256, backend version, permission, and provenance. There are currently no real Dashing scientific differences, so the status remains `conditional`.

```python
import os
from pathlib import Path

from dnakit import DNASequence
from dnakit.similarity import DashingAdapter

# The application provides the path explicitly; DNAKit does not locate or install Dashing automatically.
configured_path = os.environ.get("DNAKIT_DASHING_EXECUTABLE")
if configured_path is not None:
    dashing = DashingAdapter(Path(configured_path))
    matrix = dashing.matrix(
        (DNASequence("AACCGG"), DNASequence("AACCTT")),
        k=2,
        mode="exact",
    )
    assert matrix.provenance.backend is not None
    assert matrix.provenance.backend.name == "dashing"
```

```python
from dnakit import DNASequence
from dnakit.alignment import AlignmentConfig, align_pairwise

alignment = align_pairwise(
    DNASequence("ACGT"),
    DNASequence("ACCT"),
    config=AlignmentConfig(mode="global", gap_score=-1.0),
)
assert alignment.identity == 0.75
assert alignment.query_coverage == 1.0
```

## Dataset sorting {#datasets}

::: dnakit.datasets

Features include exact/reverse-complement/circular/IUPAC/approximate deduplication, identity/edit/k-mer/fingerprint threshold clustering, hierarchical clustering, representative sequences, random/stable hashing/hierarchical/group/similarity/temporal partitioning, joint constraint heuristic partitioning, cross-set leakage, and partition quality.

`joint_split()` does not claim global optimality; when it is not satisfied, press `infeasible_policy` to report an error or return audit results with slack marks.

## Comprehensive evaluation {#evaluation}

::: dnakit.evaluation

The reference-based method must first call `create_reference_library()` binding name, version, source, date, filter, index parameters and content digest. novelty is defined as `1 - nearest_similarity` relative to the library; memorization is exact or explicit threshold approximate replication.

The output of `evaluate_synthesis_risk()` is transparent rules and hit locations, not supplier order taking rules, structure predictions, or experiment success probabilities. `evaluate_scorecard()` Retain each component, normalization direction, weight, missing strategy and contribution.

## Molecular Biology Simulation {#molecular-biology}

::: dnakit.molbio

Digestion/end/ligation, PCR/primer matching, assembly, CRISPR, and sequence optimization are all deterministic sequence-level models. Golden Gate/BioBrick requires pre-prepared fragments with verified ends; PCR does not model yield; CRISPR does not predict efficiency or biological risk; codon optimization requires a host table from the caller.

`prepare_primer_design()` generates a validated, backend-neutral design request; the external program is executed only if the caller subsequently calls `design()` on `Primer3CLIDesignAdapter` with an explicit `primer3_core_path`. The adapter maps parameters with Boulder-IO whitelists, limits templates, number of candidates, result keys, and text length, and returns candidate coordinates, Tm, GC, penalty, warnings, permissions, and provenance. The coordinates, left and right primer sequences and product size of each candidate will be checked against the original template, and inconsistent results will be rejected based on the backend output error.

```python
from dnakit import DNASequence
from dnakit.molbio import (
    Primer3CLIDesignAdapter,
    digest_restriction,
    generate_mutation_library,
    prepare_primer_design,
)

digest = digest_restriction(DNASequence("TTTGAATTCAAA"), ("EcoRI",))
library = generate_mutation_library(
    DNASequence("ACGT"),
    {1: ("A", "G")},
    mode="single",
    seed=13,
)
assert len(digest.fragments) == 2
assert len(library.variants) == 2

request = prepare_primer_design(
    DNASequence("ACGT" * 100),
    target_start=150,
    target_end=200,
    candidate_count=2,
)
designer = Primer3CLIDesignAdapter(
    "/opt/primer3/src/primer3_core",
    thermodynamic_parameters_path="/opt/primer3/src/primer3_config",
)
if designer.info.available:
    designed = designer.design(request)
    assert designed.execution_performed
    assert designed.provenance.backend is not None
```

## Visualization {#visualization}

::: dnakit.visualization

SVG is generated natively, with all drawings using a square canvas without decorative in-figure titles; `title` in configuration is used only for SVG accessibility instructions. `save_image()` Can choose PNG, SVG or JPG, defaults to PNG when no extension is specified and `image_type` is not specified. PNG/JPG requires `viz` extra, and the original TIFF/PDF continues to be compatible. `build_html_report()` renders results explicitly provided by the caller as self-contained, filterable, expandable, read-only HTML.

## Batch processing, caching, backends and workflows {#engineering}

:::dnakit.batch

:::dnakit.cache

::: dnakit.backends

::: dnakit.workflows

Batch processing supports serial or threaded modes, stable input order, bounded in-flight, error collection, per-record seed derivation, progress callbacks, and completed ID resume. `CacheKey` Automatically include key schema and DNAKit version; caller must explicitly put input, parameters and used algorithm/backend version into components. `JSONCache` Uses normalized content keys, atomic writes, payload checksums and `max_entry_bytes` read and write caps; it has no TTL or automatic elimination.

The built-in registry has 6 stable IDs: `primer3-cli`, `nupack`, `blastn`, `mmseqs2`, `sourmash`, `dashing`. The default `primer3-cli` registration entry does not contain a path, so it remains unavailable; actual scientific calls must directly construct an explicit path adapter. BLAST/MMseqs2/sourmash's `ExternalCLIAdapter` only do passive path positioning and execute a timeout/output-capped version of the command when explicitly called by the user; they do not receive sequences, nor are they search, clustering, or sketch executors. The registry's Dashing metadata/version handle also remains passive; domain-level `DashingAdapter` is another bounded scientific computing entry that must explicitly provide an executable path and call explicitly.

The YAML/JSON pipeline uses the `dnakit-workflow-v1` strict schema and only allows `normalize/validate/descriptors/fingerprint/deduplicate/split/write/report`. Output is restricted to a dedicated directory under the configuration file; the manifest records resolved config, seed, version, step status, and artifact SHA-256. `--resume` only skips `write/report` that pass the integrity check, not a general cache. The workflow does not load arbitrary Python callables, does not execute shells/networks, and does not automatically download backends/databases.

```python
from pathlib import Path
from shutil import copy2
from tempfile import TemporaryDirectory

from dnakit import DNARecord, DNASequence
from dnakit.backends import backend_registry
from dnakit.batch import BatchConfig, run_batch
from dnakit.cache import CacheKey, JSONCache
from dnakit.workflows import load_workflow, run_workflow

tiny = [DNARecord(DNASequence("ACGT"), "a"), DNARecord(DNASequence("GGGG"), "b")]
batch = run_batch(
    tiny,
    lambda record, context: {"id": record.id, "seed": context.seed},
    name="api-example",
    config=BatchConfig(seed=7, jobs=2, execution_mode="thread", max_in_flight=2),
)
assert batch.success_count == 2
assert {
    "blastn",
    "dashing",
    "mmseqs2",
    "nupack",
    "primer3-cli",
    "sourmash",
}.issubset(set(backend_registry))

with TemporaryDirectory() as directory:
    root = Path(directory)
    cache = JSONCache(root / "cache")
    key = CacheKey.from_components(
        "example",
        {
            "sequence": "ACGT",
            "k": 2,
            "algorithm_version": "example-kmer-v1",
            "backend_version": None,
        },
    )
    cache.put(key, {"value": 1})
    assert cache.get(key) == {"value": 1}

    config_path = copy2("examples/advanced_workflow.yml", root)
    copy2("examples/fixed_demo.fasta", root)
    loaded = load_workflow(config_path)
    dry_run = run_workflow(config_path, dry_run=True)
    assert loaded.spec.schema_version == "dnakit-workflow-v1"
    assert dry_run.status == "dry-run"
```

## CLI

| Command | Function |
| ---------------------------------- | ---------------------------------- |
| `info` / `backends` | Report running environment and registered backends |
| `normalize` / `validate` | Standardization, and unified validation of a single sequence or set of records |
| `describe` | Base composition, GC, complexity and duplication reports |
| `fingerprint` | k-mer, MinHash or FracMinHash |
| `search` / `orfs` / `compare` | Pattern, ORF and sequence comparison |
| `convert` | Streaming format conversion |
| `deduplicate` / `split` | Data set organization |
| `report` | Self-contained read-only HTML report |
| `workflow` | Strict JSON/YAML multi-step workflow |

```bash
dnakit --help
dnakit workflow --help
```

`python -m dnakit.cli.workflow run CONFIG` is reserved for compatibility/development entry.

The item-by-item scope is based on the [ requirements tracking matrix ](../planning/requirements_traceability.csv).
