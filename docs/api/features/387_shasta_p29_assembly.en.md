# Shasta long-read assembly (command-line facade)

[Shasta](https://github.com/paoloshasta/shasta) assembles long Oxford Nanopore reads de novo. DNAKit runs the `shasta` command line with a fixed argument whitelist (no shell) and reads back the files Shasta writes for each pipeline stage. One function covers all stage variants of Shasta's `Assembler` API that the command line exposes: Mode 3 phased assembly from anchors, marker (k-mer) selection, MinHash alignment-candidate search, palindromic and chimeric read flagging, the alignment methods (old Shasta, SeqAn, banded SeqAn, method 4), the read graph creation methods (criteria, percentile, v3, Bayesian v4) with cross-strand edge flagging, inconsistent-alignment flagging, marker graph vertex and edge creation, duplicate-marker cleanup, transitive reduction, pruning and simplification, repeat-count consensus, detangling, and the GFA and FASTA exports.

Install Shasta separately (release binary from https://github.com/paoloshasta/shasta/releases, tested with 0.14.0) and put it on PATH as `shasta`, set `DNAKIT_SHASTA` or pass `shasta_path`. The Shasta Python API (`import shasta`) is not wrapped: it belongs to the source build, its methods only work on the in-memory state left by a previous stage, and the shared library of the release needs Python 3.10 and Ubuntu 22.04 libraries. Stages that only exist as such methods (for example `writeMarkers`, `analyzeReadGraph`, `readGraphClustering`, GFA colouring) are therefore not available; see the disposition table of the integration pack P29.

## Example data

`sim_reads.py` simulates a 150 kb random genome and long reads with 1% errors on both strands (`write_reads`). The first example uses 450 reads, the third 800 reads with 0.2% errors.

**Example file:** `sim_reads.py`

```text
import random


def write_reads(path, count=450, error=0.01, seed=1):
    """Simulate a 150 kb random genome and long reads (12-25 kb, both strands) with substitutions and indels."""
    rng = random.Random(seed)
    genome = "".join(rng.choice("ACGT") for _ in range(150_000))
    table = str.maketrans("ACGT", "TGCA")
    with open(path, "w") as out:
        for index in range(count):
            length = rng.randint(12_000, 25_000)
            start = rng.randint(0, len(genome) - length)
            read = genome[start:start + length]
            if rng.random() < 0.5:
                read = read.translate(table)[::-1]
            mutated = []
            for base in read:
                x = rng.random()
                if x < error / 3:
                    continue
                if x < 2 * error / 3:
                    mutated.append(base + rng.choice("ACGT"))
                elif x < error:
                    mutated.append(rng.choice("ACGT"))
                else:
                    mutated.append(base)
            out.write(f">r{index}\n{''.join(mutated)}\n")
```

## 1) `dnakit.assembly.shasta_assemble()`

- **Purpose:** Assemble long reads with Shasta; every stage of the pipeline is selected by a parameter (`assembly_mode`, `align_method`, `read_graph_creation_method`, `detangle_method`, ...) and any other Shasta 0.14.0 option can be given in `options`.
- **Method:** Runs `shasta --command assemble` with a whitelist of named parameters plus checked `options`; a parameter is passed only when it is given, so the chosen `config` decides otherwise. The run folder must not exist; a missing program raises `BackendUnavailableError` with an install hint.
- **API:** `dnakit.assembly.shasta_assemble(reads: 'Sequence[str]', output_directory: 'str', *, config: 'str' = 'Nanopore-May2022', assembly_mode: 'int | None' = None, anchors_path: 'str | None' = None, min_read_length: 'int | None' = None, desired_coverage: 'str | int | None' = None, palindromic_flagging: 'bool | None' = None, kmer_generation_method: 'int | None' = None, kmer_length: 'int | None' = None, kmer_probability: 'float | None' = None, minhash_m: 'int | None' = None, minhash_hash_fraction: 'float | None' = None, minhash_iteration_count: 'int | None' = None, minhash_all_pairs: 'bool | None' = None, align_method: 'int | None' = None, read_graph_creation_method: 'int | None' = None, strand_separation_method: 'int | None' = None, max_chimeric_read_distance: 'int | None' = None, flag_inconsistent_alignments: 'bool | None' = None, marker_graph_min_coverage: 'int | None' = None, marker_graph_max_coverage: 'int | None' = None, marker_graph_min_edge_coverage: 'int | None' = None, cleanup_duplicate_markers: 'bool | None' = None, simplify_max_length: 'str | Sequence[int] | None' = None, detangle_method: 'int | None' = None, consensus_caller: 'str | None' = None, write_reads_by_segment: 'bool | None' = None, store_coverage_data: 'bool | None' = None, memory_mode: 'str | None' = None, memory_backing: 'str | None' = None, options: 'Mapping[str, Any] | None' = None, threads: 'int | None' = None, shasta_path: 'str | None' = None, timeout_seconds: 'float' = 86400.0, verbose: bool = True) -> 'dict[str, Any]'`.
- **Input:** Read files (FASTA/FASTQ), a new output folder, a configuration name or file, and optional stage parameters.
- **Output:** Dictionary with `status`, `backend`, `version`, `command`, `return_code`, `assembly_directory`, `assembly_fasta`, `assembly_gfa`, `summary` (reads used, palindromic and chimeric reads, assembled segments, bases, N50), `outputs` and `log_tail`.

## 2) `dnakit.assembly.shasta_read_run_table()`

- **Purpose:** Read one of the tables Shasta writes into its run folder (read flags, MinHash and alignment histograms, read graph components, marker graph coverage histograms, reads by segment, Mode 3 tables).
- **Method:** Parses the CSV file of the run folder; numbers are converted, at most `max_rows` rows are returned.
- **API:** `dnakit.assembly.shasta_read_run_table(run_directory: 'str', table: 'str', *, max_rows: 'int' = 10000, verbose: bool = True) -> 'dict[str, Any]'`.
- **Input:** The run folder of a finished assembly and a table name.
- **Output:** Dictionary with `table`, `file`, `columns`, `row_count`, `rows`, `truncated`.

## 3) `dnakit.assembly.shasta_read_assembly_output()`

- **Purpose:** Read the assembled FASTA or the GFA graph that Shasta exported and summarise it (records, lengths, N50, GC, links).
- **Method:** Parses FASTA or GFA 1.0 (sequences or `LN:i` lengths); no program is run.
- **API:** `dnakit.assembly.shasta_read_assembly_output(run_directory: 'str', *, output: 'str' = 'fasta', file_name: 'str | None' = None, include_sequences: 'bool' = False, max_records: 'int' = 10000, verbose: bool = True) -> 'dict[str, Any]'`.
- **Input:** The run folder, the output kind (`fasta`, `gfa`, `gfa_both_strands`, `gfa_no_sequence`) or a plain `file_name` of the folder.
- **Output:** Dictionary with `records`, `record_count`, `total_length`, `longest`, `n50`, and `link_count`/`links` for GFA.

## Examples

### Example 1: Assemble reads and read back the run

Assemble simulated reads, then read the read table, the read graph components and the FASTA and GFA exports.

- **Example code:**

```python
import tempfile
from pathlib import Path
from sim_reads import write_reads
from dnakit.assembly import shasta_p29_cli_functions as sh

work = Path(tempfile.mkdtemp())
write_reads(work / "reads.fasta")
result = sh.shasta_assemble([str(work / "reads.fasta")], str(work / "run"), min_read_length=5000, threads=4, verbose=False)
summary = result["summary"]
print(result["status"], result["version"])
print(summary["reads_used"], summary["reads_flagged_chimeric"], summary["assembled_segments"])

reads = sh.shasta_read_run_table(result["assembly_directory"], "read_summary", verbose=False)
print(reads["columns"][:5], reads["row_count"])
components = sh.shasta_read_run_table(result["assembly_directory"], "read_graph_components", max_rows=2, verbose=False)
first = components["rows"][0]
print(components["row_count"], first["OrientedReadCount"], first["IsSelfComplementary"])
fasta = sh.shasta_read_assembly_output(result["assembly_directory"], verbose=False)
print(fasta["record_count"], fasta["total_length"] == summary["assembled_bases"], fasta["n50"] == summary["segment_n50"])
gfa = sh.shasta_read_assembly_output(result["assembly_directory"], output="gfa_both_strands", verbose=False)
print(gfa["record_count"], gfa["link_count"])
```

- **Example output:**

```text
completed 0.14.0
450 5 3
['Id', 'Name', 'RawLength', 'RleLength', 'RawOverRleLengthRatio'] 450
8 108 No
3 True True
6 0
```

### Example 2: Select pipeline stages with parameters

Alignment method, read graph creation method, detangling, reads by segment, Mode 2, and a refused option name.

- **Example code:**

```python
import tempfile
from pathlib import Path
from sim_reads import write_reads
from dnakit.assembly import shasta_p29_cli_functions as sh

work = Path(tempfile.mkdtemp())
write_reads(work / "reads.fasta")
runs = {
    "seqan-alignment+percentile-graph": {"align_method": 1, "read_graph_creation_method": 2},
    "bayesian-graph": {"read_graph_creation_method": 4, "strand_separation_method": 0},
    "detangle+reads-by-segment": {"detangle_method": 2, "write_reads_by_segment": True},
    "phased-mode-2": {"assembly_mode": 2, "strand_separation_method": 2},
}
for name, parameters in runs.items():
    result = sh.shasta_assemble([str(work / "reads.fasta")], str(work / name), min_read_length=5000, threads=4, verbose=False, **parameters)
    print(name, result["status"])
by_segment = sh.shasta_read_run_table(str(work / "detangle+reads-by-segment"), "reads_by_segment", max_rows=1, verbose=False)
print(by_segment["row_count"] > 0, by_segment["columns"][:4])
try:
    sh.shasta_assemble([str(work / "reads.fasta")], str(work / "bad"), options={"Align.noSuchOption": 1}, verbose=False)
except Exception as error:
    print(type(error).__name__)
```

- **Example output:**

```text
seqan-alignment+percentile-graph completed
bayesian-graph completed
detangle+reads-by-segment completed
phased-mode-2 completed
True ['AssembledSegmentId', 'EdgeCount', 'OrientedReadCount', 'OrientedReadId']
ConfigurationError
```

### Example 3: Mode 3 phased assembly from anchors

Mode 3 assembly with the raw-sequence configuration `Nanopore-ncm23-May2024` and the Mode 3 tables.

- **Example code:**

```python
import tempfile
from pathlib import Path
from sim_reads import write_reads
from dnakit.assembly import shasta_p29_cli_functions as sh

work = Path(tempfile.mkdtemp())
write_reads(work / "reads.fasta", count=800, error=0.002)
result = sh.shasta_assemble(
    [str(work / "reads.fasta")], str(work / "mode3"), config="Nanopore-ncm23-May2024", assembly_mode=3,
    min_read_length=5000, threads=4, verbose=False,
)
print(result["return_code"], result["version"])
print(sorted(Path(result["assembly_directory"]).glob("Anchor*"))[0].name)
anchors = sh.shasta_read_run_table(result["assembly_directory"], "anchor_coverage_histogram", verbose=False)
print(anchors["columns"], anchors["row_count"] > 0)
print(sh.shasta_read_run_table(result["assembly_directory"], "components", max_rows=1, verbose=False)["columns"][:3])
```

- **Example output:**

```text
0 0.14.0
AnchorCoverageHistogram.csv
['Coverage', 'Frequency'] True
['Component', 'Reads', 'Segments']
```
