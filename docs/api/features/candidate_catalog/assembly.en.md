# Assembly

This category contains **37** documented capabilities. 27 compute a result natively in pure Python from the supplied input (the 5 contiguity statistics `calculate_l50`, `calculate_n50`, `calculate_ng50`, `calculate_nx`, `calculate_stats`, whose sections are unchanged, and 22 functions for pangenome matrices, k-mer models, assembly checks and comparisons, repeats and telomeres, GFA graphs, parsing, scaffolding with an explicit layout and contaminant removal). 9 always return `status: "not_implemented"` because they need reads, alignments, a reference or a database that is not part of the input (`build_chromosome_level_assembly`, `build_haplotype_assembly`, `detect_assembly_contamination`, `detect_chimeric_contigs`, `detect_misassemblies`, `identify_centromeres`, `phase_assembly`, `polish`, `purge_haplotigs`). `run` executes a caller-installed assembler (Flye, hifiasm, MEGAHIT, Shasta or Unicycler) and returns `not_implemented` while that program is not on PATH. The documented functions are exported through the formal DNAKit namespace; functions with a compatible JSON call plan are also exposed through the Agent/MCP registry.

## 1) add_genome_to_pangenome

- **Purpose:** Adds one genome's gene set as a new column of a gene presence/absence matrix.

- **Method:** Reads the existing pangenome matrix, appends a column for the new genome (1 where the gene is in `genes`, else 0) and appends an all-absent row for the old genomes for every gene not seen before (new genes are sorted alphabetically). No homology search is done: a gene name is its own cluster, so genes that are the same ortholog under different names count as new. The new genome name must not already exist.

- **API:** `dnakit.assembly.add_genome_to_pangenome(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input` (required): gene presence/absence data: `{genome: [genes]}`, `{gene: {genome: 0/1}}` or `{"genomes": [...], "genes": [...], "matrix": [[...]]}` (rows = genes, columns = genomes; values are booleans or non-negative numbers, > 0 meaning present). Gene names are taken as already-clustered orthology groups; no sequence clustering or homology search is done. `options` (required keys): `genome` (str, the new unique genome name), `genes` (list of str, the genes present in the new genome).

- **Output:** Standard envelope: `status` (`"ok"`), `module`, `operation`, `provenance`, `input_summary`, `warnings` and `result`. `result` contains `genomes`, `genes`, `matrix` (rows = genes, columns = genomes, values 0/1) of the updated pangenome, `new_gene_count`, `new_genes`, `pangenome_size_before` and `pangenome_size_after`.

- **Example code:**

```python
import json

from dnakit.assembly import add_genome_to_pangenome

pangenome = {"A": ["g1", "g2"], "B": ["g2", "g3"]}
result = add_genome_to_pangenome(pangenome, options={"genome": "C", "genes": ["g2", "g4"]})
print(json.dumps(result, indent=2))
```

- **Example output:**

```json
{
  "module": "assembly",
  "operation": "add_genome_to_pangenome",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 1,
    "text_length": 38,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "6d290492cc572135bfe94b6f4b9597f82691599f0ace9fd6682f9798aab95803"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "genomes": ["A", "B", "C"],
    "genes": ["g1", "g2", "g3", "g4"],
    "matrix": [
      [1, 0, 0],
      [1, 1, 1],
      [0, 1, 0],
      [0, 0, 1]
    ],
    "new_gene_count": 1,
    "new_genes": ["g4"],
    "pangenome_size_before": 3,
    "pangenome_size_after": 4
  }
}
```

## 2) build_chromosome_level_assembly

- **Purpose:** Chromosome-level assembly is not implemented: ordering and orienting contigs into chromosomes needs Hi-C, optical-map or synteny evidence.

- **Method:** Not implemented. Joining contigs into chromosomes needs linking evidence (Hi-C contacts, optical maps or a reference genome) and a scaffolder such as YaHS or SALSA, none of which is part of the input; use `scaffold` with an explicit layout if the order and orientation are already known. No heuristic stand-in is applied, so the function never fabricates a result; it validates nothing about the input beyond summarising it in `input_summary` and returns a structured `not_implemented` envelope.

- **API:** `dnakit.assembly.build_chromosome_level_assembly(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input` (required): anything (it is only summarised in `input_summary`). `options` (optional): ignored.

- **Output:** Always returns `status: "not_implemented"` (with `module`, `operation`, `provenance`, `input_summary`, `warnings`), `required_dependency` (what is missing: scaffolding evidence (Hi-C contacts, optical maps or a reference) and a scaffolder such as YaHS/SALSA) and `message` (why it cannot be computed from the input). No `result` is produced.

- **Example code:**

```python
import json

from dnakit.assembly import build_chromosome_level_assembly

result = build_chromosome_level_assembly(">c1\nACGTACGT\n>c2\nGGCCGGCC\n")
print(json.dumps(result, indent=2))
```

- **Example output:**

```json
{
  "module": "assembly",
  "operation": "build_chromosome_level_assembly",
  "provenance": {
    "implementation": "not_implemented",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 2,
    "text_length": 26,
    "sequence_bases": 16,
    "lengths": [8, 8],
    "sha256": "eaa8633c7014be534dfd6571d98374f8f6d24defb95b394e7f9ad5bbe9ac56f4"
  },
  "warnings": [],
  "status": "not_implemented",
  "required_dependency": "scaffolding evidence (Hi-C contacts, optical maps or a reference) and a scaffolder such as YaHS/SALSA",
  "message": "Ordering and orienting contigs into chromosomes needs linking evidence that is not part of the input; use scaffold() with an explicit layout instead."
}
```

## 3) build_command

- **Purpose:** Builds, without running it, the command line of one supported assembler (Flye, hifiasm, MEGAHIT, Shasta or Unicycler) as an argument list.

- **Method:** Assembles the argument vector with the same rules as the `assemble_flye`, `assemble_hifiasm`, `assemble_megahit`, `assemble_shasta` and `assemble_unicycler` wrappers, validating tool-specific required keys (for example Flye needs `reads` and `genome_size`; Unicycler needs `read1` or `long_reads`) and value ranges. Nothing is executed, no file is read or created and the program does not have to be installed.

- **API:** `dnakit.assembly.build_command(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input` (required): a mapping of settings; `options` are merged over it and win. Keys: `tool` (required: `flye`, `hifiasm`, `megahit`, `shasta` or `unicycler`), `output_directory` (required), `program` (executable name, default = the tool name), `threads` (1-256, default 1), `reads` (flye, hifiasm, shasta), `read1` / `read2` (megahit, unicycler), `long_reads` (unicycler), `genome_size` (flye, for example `"5m"`), `platform` (flye, default `"nano-hq"`), `config` (shasta, default `"Nanopore-May2022"`), `mode` (unicycler, default `"normal"`), `prefix` (hifiasm, default `"assembly"`), `min_contig_length` (megahit).

- **Output:** Standard envelope: `status` (`"ok"`), `module`, `operation`, `provenance`, `input_summary`, `warnings` and `result`. `result` contains `tool`, `argv` (the argument list, program name first) and `command_line` (the shell-quoted command string).

- **Example code:**

```python
import json

from dnakit.assembly import build_command

settings = {
    "tool": "flye",
    "reads": "reads.fastq",
    "genome_size": "5m",
    "output_directory": "flye_out",
    "threads": 4,
}
result = build_command(settings)
print(json.dumps(result, indent=2))
```

- **Example output:**

```json
{
  "module": "assembly",
  "operation": "build_command",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 1,
    "text_length": 107,
    "sequence_bases": 10,
    "lengths": [10],
    "sha256": "b9a386c642fbd16778c7d5400bc52c7bba3407abbe68afa7505c06e2602b5108"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "tool": "flye",
    "argv": ["flye", "--nano-hq", "reads.fastq", "--genome-size", "5m", "--out-dir", "flye_out", "--threads", "4"],
    "command_line": "flye --nano-hq reads.fastq --genome-size 5m --out-dir flye_out --threads 4"
  }
}
```

## 4) build_graph_genome

- **Purpose:** Builds a compacted forward-strand de Bruijn graph (unitigs, written as GFA1) from input sequences.

- **Method:** Every k-mer of every sequence (split at characters other than A/C/G/T) is a node; an edge joins two k-mers observed consecutively, i.e. with a (k-1)-base overlap. Unitigs are the maximal non-branching paths of this graph (isolated cycles are also reported); each unitig records the input sequences it occurs in (`sources`). The reverse-complement strand is NOT merged (the graph is not bidirected, so a sequence and its reverse complement give separate nodes), and no alignment or variation-graph construction is performed. Input is limited to 2,000,000 bases.

- **API:** `dnakit.assembly.build_graph_genome(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input` (required): FASTA text or path, `{name: sequence}` or a list of sequences (each is one input sequence, named `contig_<i>` when unnamed). `options` (optional): `k` (int 3-127, default 21), `output` (optional file path to write the GFA1 text).

- **Output:** Standard envelope: `status` (`"ok"`), `module`, `operation`, `provenance`, `input_summary`, `warnings` and `result`. `result` contains `k`, `kmer_count`, `unitig_count`, `link_count`, `unitigs` (each with `name` `utg<i>`, `sequence`, `length`, `kmers`, `sources`), `links` (`from`, `to`, `overlap`), `gfa` (GFA1 text: `S` lines with an `LN` tag and `L` lines with a `(k-1)M` overlap, all `+`), `note` and `output_path`.

- **Example code:**

```python
import json

from dnakit.assembly import build_graph_genome

result = build_graph_genome({"s1": "ACGTT", "s2": "ACGTA"}, options={"k": 3})
print(json.dumps(result, indent=2))
```

- **Example output:**

```json
{
  "module": "assembly",
  "operation": "build_graph_genome",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 1,
    "text_length": 30,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "db2285ea0f9ff1176712fd147fdb9ef5c1c1bf968cc8ca6a4a4ecdbe64e768ad"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "k": 3,
    "kmer_count": 4,
    "unitig_count": 3,
    "link_count": 2,
    "unitigs": [
      {
        "name": "utg1",
        "sequence": "ACGT",
        "length": 4,
        "kmers": 2,
        "sources": ["s1", "s2"]
      },
      {
        "name": "utg2",
        "sequence": "GTA",
        "length": 3,
        "kmers": 1,
        "sources": ["s2"]
      },
      {
        "name": "utg3",
        "sequence": "GTT",
        "length": 3,
        "kmers": 1,
        "sources": ["s1"]
      }
    ],
    "links": [
      {
        "from": "utg1",
        "to": "utg2",
        "overlap": 2
      },
      {
        "from": "utg1",
        "to": "utg3",
        "overlap": 2
      }
    ],
    "gfa": "H\tVN:Z:1.0\nS\tutg1\tACGT\tLN:i:4\nS\tutg2\tGTA\tLN:i:3\nS\tutg3\tGTT\tLN:i:3\nL\tutg1\t+\tutg2\t+\t2M\nL\tutg1\t+\tutg3\t+\t2M\n",
    "note": "forward-strand (non-canonical) de Bruijn graph; nodes are k-mers, edges are observed (k-1)-overlaps; unitigs are maximal non-branching paths",
    "output_path": null
  }
}
```

## 5) build_haplotype_assembly

- **Purpose:** Haplotype-resolved assembly is not implemented: it needs reads (HiFi/ONT, Hi-C or parental k-mers) and a phasing assembler.

- **Method:** Not implemented. A haplotype-resolved assembly cannot be derived from contigs alone; it needs sequencing reads and an assembler with phasing support such as hifiasm (see `build_command` / `run`). No heuristic stand-in is applied, so the function never fabricates a result; it validates nothing about the input beyond summarising it in `input_summary` and returns a structured `not_implemented` envelope.

- **API:** `dnakit.assembly.build_haplotype_assembly(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input` (required): anything (it is only summarised in `input_summary`). `options` (optional): ignored.

- **Output:** Always returns `status: "not_implemented"` (with `module`, `operation`, `provenance`, `input_summary`, `warnings`), `required_dependency` (what is missing: an assembler with phasing support (for example hifiasm) and sequencing reads) and `message` (why it cannot be computed from the input). No `result` is produced.

- **Example code:**

```python
import json

from dnakit.assembly import build_haplotype_assembly

result = build_haplotype_assembly(">c1\nACGTACGT\n>c2\nGGCCGGCC\n")
print(json.dumps(result, indent=2))
```

- **Example output:**

```json
{
  "module": "assembly",
  "operation": "build_haplotype_assembly",
  "provenance": {
    "implementation": "not_implemented",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 2,
    "text_length": 26,
    "sequence_bases": 16,
    "lengths": [8, 8],
    "sha256": "eaa8633c7014be534dfd6571d98374f8f6d24defb95b394e7f9ad5bbe9ac56f4"
  },
  "warnings": [],
  "status": "not_implemented",
  "required_dependency": "an assembler with phasing support (for example hifiasm) and sequencing reads",
  "message": "Haplotype-resolved assembly cannot be derived from contigs alone; run assemble_hifiasm on reads instead."
}
```

## 6) build_pangenome

- **Purpose:** Builds a gene presence/absence pangenome matrix with core, shell and cloud counts from per-genome gene sets.

- **Method:** Collects the union of genes into a genes x genomes 0/1 matrix and partitions the genes by the fraction of genomes that carry them: core = present in at least `ceil(core_threshold * genomes)` genomes; of the remaining accessory genes, cloud = present in a fraction below `cloud_threshold`, shell = the rest. Also counts genes and genome-unique genes per genome. Gene names are taken as already-clustered orthology groups; no sequence clustering is done.

- **API:** `dnakit.assembly.build_pangenome(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input` (required): gene presence/absence data: `{genome: [genes]}`, `{gene: {genome: 0/1}}` or `{"genomes": [...], "genes": [...], "matrix": [[...]]}` (rows = genes, columns = genomes; values are booleans or non-negative numbers, > 0 meaning present). Gene names are taken as already-clustered orthology groups; no sequence clustering or homology search is done. `options` (optional): `core_threshold` (fraction in (0, 1], default 0.99), `cloud_threshold` (fraction in [0, 1], default 0.15).

- **Output:** Standard envelope: `status` (`"ok"`), `module`, `operation`, `provenance`, `input_summary`, `warnings` and `result`. `result` contains `genomes`, `genes`, `matrix`, `pangenome_size`, `core_count`, `shell_count`, `cloud_count`, `core_threshold`, `genes_per_genome` and `unique_genes_per_genome`.

- **Example code:**

```python
import json

from dnakit.assembly import build_pangenome

genomes = {"G1": ["a", "b", "c"], "G2": ["a", "b", "d"], "G3": ["a", "e"]}
result = build_pangenome(genomes, options={"cloud_threshold": 0.5})
print(json.dumps(result, indent=2))
```

- **Example output:**

```json
{
  "module": "assembly",
  "operation": "build_pangenome",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 1,
    "text_length": 64,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "e4372679c8d0e6836c36cd361cef3caf9656e4da871974b5c7cf4f7878ed448b"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "genomes": ["G1", "G2", "G3"],
    "genes": ["a", "b", "c", "d", "e"],
    "matrix": [
      [1, 1, 1],
      [1, 1, 0],
      [1, 0, 0],
      [0, 1, 0],
      [0, 0, 1]
    ],
    "pangenome_size": 5,
    "core_count": 1,
    "shell_count": 1,
    "cloud_count": 3,
    "core_threshold": 0.99,
    "genes_per_genome": {
      "G1": 3,
      "G2": 3,
      "G3": 2
    },
    "unique_genes_per_genome": {
      "G1": 1,
      "G2": 1,
      "G3": 1
    }
  }
}
```

## 7) calculate_assembly_completeness

- **Purpose:** Estimates assembly completeness as a length ratio and/or as the recovery of the reference's k-mers.

- **Method:** Length mode: `length_ratio = assembly length / reference_length` and `length_completeness = min(ratio, 1)`. K-mer mode: `kmer_completeness` = the fraction of the reference's distinct canonical k-mers (a k-mer and its reverse complement count as one; only A/C/G/T stretches are used) that also occur in the assembly. This is a length / k-mer measure, not a gene-content measure (see `calculate_busco_completeness`). At least one of `reference_length` or `reference` is required; k-mer mode needs sequences and is limited to 5,000,000 bases in total.

- **API:** `dnakit.assembly.calculate_assembly_completeness(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input` (required): an assembly given as FASTA text or a path to a FASTA file, `{name: sequence}`, or a list of sequences / `{"name", "sequence"}` records (contig lengths, as `{name: length}` or a list of integers, are also accepted). `options` (at least one of `reference_length` / `reference`): `reference_length` (positive number, expected genome size), `reference` (an assembly with sequences, in any of the forms above), `k` (int 3-63, default 21; k-mer mode).

- **Output:** Standard envelope: `status` (`"ok"`), `module`, `operation`, `provenance`, `input_summary`, `warnings` and `result`. `result` contains `assembly_length`, `contig_count`, `methods` (`"length"` and/or `"kmer"`) and, per mode, `reference_length`, `length_ratio`, `length_completeness`; `k`, `reference_distinct_kmers`, `shared_kmers`, `kmer_completeness`.

- **Example code:**

```python
import json

from dnakit.assembly import calculate_assembly_completeness

result = calculate_assembly_completeness(
    {"c1": "ACGTA"},
    options={"reference_length": 8, "reference": "ACGTACGT", "k": 4},
)
print(json.dumps(result, indent=2))
```

- **Example output:**

```json
{
  "module": "assembly",
  "operation": "calculate_assembly_completeness",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 1,
    "text_length": 15,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "6fbe1ad637a90ef6b20a533c11a2dfc6669bbeaa5ffcabfcfda1d2880ded589e"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "assembly_length": 5,
    "contig_count": 1,
    "reference_length": 8,
    "length_ratio": 0.625,
    "length_completeness": 0.625,
    "k": 4,
    "reference_distinct_kmers": 3,
    "shared_kmers": 2,
    "kmer_completeness": 0.6666666666666666,
    "methods": ["length", "kmer"]
  }
}
```

## 8) calculate_busco_completeness

- **Purpose:** Computes BUSCO completeness percentages from BUSCO counts or from a BUSCO short-summary line (it does not run BUSCO).

- **Method:** complete = single (S) + duplicated (D); each percentage is the count divided by n, where n is the stated total (summary text) or S + D + F + M (counts). From a summary line such as `C:98.5%[S:97.0%,D:1.5%],F:0.5%,M:1.0%,n:1614` the counts are recovered as `round(percent * n / 100)`. It only summarises existing BUSCO results; to run BUSCO itself use `run_busco`.

- **API:** `dnakit.assembly.calculate_busco_completeness(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input` (required): a mapping with integer counts for single (`single`, or `S`), duplicated (`duplicated` / `D`), fragmented (`fragmented` / `F`) and missing (`missing` / `M`), or the text or path of a BUSCO short summary containing the `C:..%[S:..%,D:..%],F:..%,M:..%,n:..` line. `options`: none.

- **Output:** Standard envelope: `status` (`"ok"`), `module`, `operation`, `provenance`, `input_summary`, `warnings` and `result`. `result` contains `n`, `single`, `duplicated`, `complete`, `fragmented`, `missing`, `complete_percent`, `single_percent`, `duplicated_percent`, `fragmented_percent`, `missing_percent` and `summary` (the BUSCO-style summary string).

- **Example code:**

```python
import json

from dnakit.assembly import calculate_busco_completeness

summary = "# BUSCO\n\tC:98.5%[S:97.0%,D:1.5%],F:0.5%,M:1.0%,n:1614\n"
result = calculate_busco_completeness(summary)
print(json.dumps(result, indent=2))
```

- **Example output:**

```json
{
  "module": "assembly",
  "operation": "calculate_busco_completeness",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 2,
    "text_length": 54,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "b543a46533d7392dd244a0fae5395eb128a9c9e83904ce3338824e9c42dda654"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "n": 1614,
    "single": 1566,
    "duplicated": 24,
    "complete": 1590,
    "fragmented": 8,
    "missing": 16,
    "complete_percent": 98.51301115241635,
    "single_percent": 97.02602230483271,
    "duplicated_percent": 1.486988847583643,
    "fragmented_percent": 0.49566294919454773,
    "missing_percent": 0.9913258983890955,
    "summary": "C:98.5%[S:97.0%,D:1.5%],F:0.5%,M:1.0%,n:1614"
  }
}
```

## 9) calculate_l50

- **Purpose:** Document a reusable calculate l50 capability for Assembly workflows. Candidate CAP-0061 comes from 1 source row(s) and is currently marked `new_candidate`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the calculate l50 operation, and return structured results with provenance and diagnostics. The dependency-free implementation is available at `src/dnakit/assembly/local.py` and is exported through the formal DNAKit namespace. It computes deterministic contiguity statistics from explicit contig lengths; no assembler, graph, reference genome, or BUSCO database is queried. Review basis: `no_exact_public_match`. No exact public-function name was found in the comparison.

- **API:** `dnakit.assembly.calculate_l50(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `input`; optional arguments: `options`. IO arguments accept text, bytes, mappings, paths, or materialized records as required by the signature.

- **Output:** The local `calculate_l50` function returns a completed JSON-serializable envelope with `status`, `operation`, `result`, and explicit local provenance. The metric is computed from positive contig lengths only; no assembler, graph, reference genome, or BUSCO database is queried.

- **Example code:**

```python
from pathlib import Path

from dnakit.assembly import calculate_l50

result = calculate_l50(
    input={"lengths": [100, 60, 40]},
    options={}
)
print(result)
```

- **Example output:**

```json
{
  "status": "completed",
  "capability_id": "CAP-0061",
  "operation": "calculate_l50",
  "result_type": "dict",
  "result_preview": "{'status': 'completed', 'operation': 'calculate_l50', 'result': {'l50': 1, 'n_at_percentile': 100, 'percentile': 50.0, 'contig_count': 3, 'total_length': 200}, 'provenance': {'implementation': 'dnakit.native_assembly_metrics', 'backend': 'local_dependency_free', 'scientific_boundary': 'contiguity...",
  "provenance": {
    "implementation": "dnakit.native_assembly_metrics",
    "backend": "local_dependency_free",
    "scientific_boundary": "contiguity statistics over explicit contig lengths; no assembler or reference genome"
  }
}
```

## 10) calculate_n50

- **Purpose:** Document a reusable calculate n50 capability for Assembly workflows. Candidate CAP-0062 comes from 1 source row(s) and is currently marked `new_candidate`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the calculate n50 operation, and return structured results with provenance and diagnostics. The dependency-free implementation is available at `src/dnakit/assembly/local.py` and is exported through the formal DNAKit namespace. It computes deterministic contiguity statistics from explicit contig lengths; no assembler, graph, reference genome, or BUSCO database is queried. Review basis: `no_exact_public_match`. No exact public-function name was found in the comparison.

- **API:** `dnakit.assembly.calculate_n50(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `input`; optional arguments: `options`. IO arguments accept text, bytes, mappings, paths, or materialized records as required by the signature.

- **Output:** The local `calculate_n50` function returns a completed JSON-serializable envelope with `status`, `operation`, `result`, and explicit local provenance. The metric is computed from positive contig lengths only; no assembler, graph, reference genome, or BUSCO database is queried.

- **Example code:**

```python
from pathlib import Path

from dnakit.assembly import calculate_n50

result = calculate_n50(
    input={"lengths": [100, 60, 40]},
    options={}
)
print(result)
```

- **Example output:**

```json
{
  "status": "completed",
  "capability_id": "CAP-0062",
  "operation": "calculate_n50",
  "result_type": "dict",
  "result_preview": "{'status': 'completed', 'operation': 'calculate_n50', 'result': {'n50': 100, 'l50': 1, 'contig_count': 3, 'total_length': 200, 'percentile': 50.0, 'options_used': []}, 'provenance': {'implementation': 'dnakit.native_assembly_metrics', 'backend': 'local_dependency_free', 'scientific_boundary': 'co...",
  "provenance": {
    "implementation": "dnakit.native_assembly_metrics",
    "backend": "local_dependency_free",
    "scientific_boundary": "contiguity statistics over explicit contig lengths; no assembler or reference genome"
  }
}
```

## 11) calculate_ng50

- **Purpose:** Document a reusable calculate ng50 capability for Assembly workflows. Candidate CAP-0063 comes from 1 source row(s) and is currently marked `new_candidate`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the calculate ng50 operation, and return structured results with provenance and diagnostics. The dependency-free implementation is available at `src/dnakit/assembly/local.py` and is exported through the formal DNAKit namespace. It computes deterministic contiguity statistics from explicit contig lengths; no assembler, graph, reference genome, or BUSCO database is queried. Review basis: `no_exact_public_match`. No exact public-function name was found in the comparison.

- **API:** `dnakit.assembly.calculate_ng50(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `input`; optional arguments: `options`. IO arguments accept text, bytes, mappings, paths, or materialized records as required by the signature.

- **Output:** The local `calculate_ng50` function returns a completed JSON-serializable envelope with `status`, `operation`, `result`, and explicit local provenance. The metric is computed from positive contig lengths only; no assembler, graph, reference genome, or BUSCO database is queried.

- **Example code:**

```python
from pathlib import Path

from dnakit.assembly import calculate_ng50

result = calculate_ng50(
    input={"lengths": [100, 60, 40]},
    options={"genome_size": 250}
)
print(result)
```

- **Example output:**

```json
{
  "status": "completed",
  "capability_id": "CAP-0063",
  "operation": "calculate_ng50",
  "result_type": "dict",
  "result_preview": "{'status': 'completed', 'operation': 'calculate_ng50', 'result': {'ng': 60, 'ng50': 60, 'contigs_to_threshold': 2, 'percentile': 50.0, 'expected_genome_size': 250, 'assembly_length': 200}, 'provenance': {'implementation': 'dnakit.native_assembly_metrics', 'backend': 'local_dependency_free', 'scie...",
  "provenance": {
    "implementation": "dnakit.native_assembly_metrics",
    "backend": "local_dependency_free",
    "scientific_boundary": "contiguity statistics over explicit contig lengths; no assembler or reference genome"
  }
}
```

## 12) calculate_nx

- **Purpose:** Document a reusable calculate nx capability for Assembly workflows. Candidate CAP-0064 comes from 1 source row(s) and is currently marked `new_candidate`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the calculate nx operation, and return structured results with provenance and diagnostics. The dependency-free implementation is available at `src/dnakit/assembly/local.py` and is exported through the formal DNAKit namespace. It computes deterministic contiguity statistics from explicit contig lengths; no assembler, graph, reference genome, or BUSCO database is queried. Review basis: `no_exact_public_match`. Existing-name evidence: `calculate_dn_ds`.

- **API:** `dnakit.assembly.calculate_nx(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `input`; optional arguments: `options`. IO arguments accept text, bytes, mappings, paths, or materialized records as required by the signature.

- **Output:** The local `calculate_nx` function returns a completed JSON-serializable envelope with `status`, `operation`, `result`, and explicit local provenance. The metric is computed from positive contig lengths only; no assembler, graph, reference genome, or BUSCO database is queried.

- **Example code:**

```python
from pathlib import Path

from dnakit.assembly import calculate_nx

result = calculate_nx(
    input={"lengths": [100, 60, 40]},
    options={}
)
print(result)
```

- **Example output:**

```json
{
  "status": "completed",
  "capability_id": "CAP-0064",
  "operation": "calculate_nx",
  "result_type": "dict",
  "result_preview": "{'status': 'completed', 'operation': 'calculate_nx', 'result': {'percentile': 50.0, 'nx': 100, 'contigs_to_threshold': 1, 'contig_count': 3, 'total_length': 200}, 'provenance': {'implementation': 'dnakit.native_assembly_metrics', 'backend': 'local_dependency_free', 'scientific_boundary': 'contigu...",
  "provenance": {
    "implementation": "dnakit.native_assembly_metrics",
    "backend": "local_dependency_free",
    "scientific_boundary": "contiguity statistics over explicit contig lengths; no assembler or reference genome"
  }
}
```

## 13) calculate_pangenome_openness

- **Purpose:** Fits Heaps' law to the number of new genes per added genome to classify a pangenome as open or closed.

- **Method:** Model (Tettelin et al. 2008): the mean number of new genes found when the N-th genome is added is `kappa * N^(-alpha)`. The genomes are added in `n_permutations` random orders (a fixed `seed` makes it deterministic; `0` uses the given order once), the new-gene counts are averaged per N, and `alpha`, `kappa` are fitted by ordinary least squares on log-log for N = 2..G (N values with zero new genes are dropped). The pangenome is open if `alpha < 1` and closed if `alpha >= 1`; the `model` string in the result states the same rule. At least four genomes are needed, and at least three N values with new genes. The first genome (N = 1) is not part of the fit.

- **API:** `dnakit.assembly.calculate_pangenome_openness(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input` (required): gene presence/absence data: `{genome: [genes]}`, `{gene: {genome: 0/1}}` or `{"genomes": [...], "genes": [...], "matrix": [[...]]}` (rows = genes, columns = genomes; values are booleans or non-negative numbers, > 0 meaning present). Gene names are taken as already-clustered orthology groups; no sequence clustering or homology search is done (at least 4 genomes). `options` (optional): `n_permutations` (int 0-10,000, default 100), `seed` (int, default 0).

- **Output:** Standard envelope: `status` (`"ok"`), `module`, `operation`, `provenance`, `input_summary`, `warnings` and `result`. `result` contains `alpha`, `kappa`, `points_used`, `points_dropped`, `open` (true if `alpha < 1`), `genome_count`, `permutations`, `seed` (null when `n_permutations` is 0), `mean_new_genes_by_genome_number`, `mean_pangenome_size_by_genome_number` and a `model` string.

- **Example code:**

```python
import json

from dnakit.assembly import calculate_pangenome_openness

genomes = {
    "g1": list("abcdefgh"),
    "g2": list("abcd") + list("ijkl"),
    "g3": list("ab") + list("mn"),
    "g4": ["a", "o"],
    "g5": ["a", "p"],
}
result = calculate_pangenome_openness(genomes, options={"n_permutations": 0})
print(json.dumps(result, indent=2))
```

- **Example output:**

```json
{
  "module": "assembly",
  "operation": "calculate_pangenome_openness",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 1,
    "text_length": 160,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "e6415dfbb17e1fa60f462f0dea884a8f08171b89837cce9dfad1af2832d0e517"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "alpha": 1.6324718703854002,
    "kappa": 11.866529443446138,
    "points_used": 4,
    "points_dropped": 0,
    "open": false,
    "genome_count": 5,
    "permutations": 1,
    "seed": null,
    "mean_new_genes_by_genome_number": [8.0, 4.0, 2.0, 1.0, 1.0],
    "mean_pangenome_size_by_genome_number": [8.0, 12.0, 14.0, 15.0, 16.0],
    "model": "Heaps' law: mean number of new genes found when the N-th genome is added = kappa * N^(-alpha), fitted by OLS on log-log for N=2..G; open if alpha < 1, closed if alpha >= 1 (Tettelin et al. 2008; Rasko et al. 2008)"
  }
}
```

## 14) calculate_stats

- **Purpose:** Document a reusable calculate stats capability for Assembly workflows. Candidate CAP-0066 comes from 10 source row(s) and is currently marked `review_needed`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the calculate stats operation, and return structured results with provenance and diagnostics. The optional `assembler` parameter selects one supported variant (the supported values recorded in the source catalog). The dependency-free implementation is available at `src/dnakit/assembly/local.py` and is exported through the formal DNAKit namespace. It computes deterministic contiguity statistics from explicit contig lengths; no assembler, graph, reference genome, or BUSCO database is queried. Review basis: `merged_dimension`. No exact public-function name was found in the comparison.

- **API:** `dnakit.assembly.calculate_stats(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `input`; optional arguments: `options`. IO arguments accept text, bytes, mappings, paths, or materialized records as required by the signature.

- **Output:** The local `calculate_stats` function returns a completed JSON-serializable envelope with `status`, `operation`, `result`, and explicit local provenance. The metric is computed from positive contig lengths only; no assembler, graph, reference genome, or BUSCO database is queried.

- **Example code:**

```python
from pathlib import Path

from dnakit.assembly import calculate_stats

result = calculate_stats(
    input={"lengths": [100, 60, 40]},
    options={}
)
print(result)
```

- **Example output:**

```json
{
  "status": "completed",
  "capability_id": "CAP-0066",
  "operation": "calculate_stats",
  "result_type": "dict",
  "result_preview": "{'status': 'completed', 'operation': 'calculate_stats', 'result': {'contig_count': 3, 'total_length': 200, 'min_length': 40, 'max_length': 100, 'mean_length': 66.66666666666667, 'n50': 100, 'l50': 1, 'n90': 40, 'l90': 3, 'sorted_lengths': [100, 60, 40], 'options_used': []}, 'provenance': {'implem...",
  "provenance": {
    "implementation": "dnakit.native_assembly_metrics",
    "backend": "local_dependency_free",
    "scientific_boundary": "contiguity statistics over explicit contig lengths; no assembler or reference genome"
  }
}
```

## 15) check_input

- **Purpose:** Validates assembly input and reports the detected format, contig count, total length, duplicate contig names and non-IUPAC symbols.

- **Method:** Parses the input as an assembly (sequences or only lengths), then checks that contig names are unique and that every sequence consists only of the IUPAC nucleotide symbols `ACGTNRYSWKMBDHV` (case-insensitive; `-` and `*` are reported as invalid). Unusable input (empty, wrong type, not FASTA) is reported as `valid: false` with the reason in `errors` instead of raising.

- **API:** `dnakit.assembly.check_input(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input` (required): an assembly given as FASTA text or a path to a FASTA file, `{name: sequence}`, or a list of sequences / `{"name", "sequence"}` records (contig lengths, as `{name: length}` or a list of integers, are also accepted). `options`: none.

- **Output:** Standard envelope: `status` (`"ok"`), `module`, `operation`, `provenance`, `input_summary`, `warnings` and `result`. `result` contains `valid`, `format` (`"sequences"`, `"lengths"`, or null when unusable), `contig_count`, `total_length`, `duplicate_names`, `invalid_symbols` (contig -> offending symbols) and `errors` (list of messages); for unusable input only `valid`, `format`, `contig_count` and `errors`.

- **Example code:**

```python
import json

from dnakit.assembly import check_input

result = check_input(">c1\nACGT\n>c1\nAXGT\n")
print(json.dumps(result, indent=2))
```

- **Example output:**

```json
{
  "module": "assembly",
  "operation": "check_input",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 2,
    "text_length": 18,
    "sequence_bases": 8,
    "lengths": [4, 4],
    "sha256": "0681afaaa49db2047869730361f15d87a0c5b07c6027752a89404eba33894188"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "valid": false,
    "format": "sequences",
    "contig_count": 2,
    "total_length": 8,
    "duplicate_names": ["c1"],
    "invalid_symbols": {
      "c1": ["X"]
    },
    "errors": ["duplicate contig names: c1", "non-IUPAC symbols in 1 contig(s) (c1: X)"]
  }
}
```

## 16) compare_assemblies

- **Purpose:** Compares two assemblies by contiguity statistics (count, length, N50, auN, GC, gaps) and by the contig sequences they share.

- **Method:** For each assembly it computes the contig count, total/min/max/mean length, N50/L50 and N90/L90 (longest-first cumulative length reaching 50% / 90% of the assembly), auN = sum(L_i^2) / sum(L_i) (area under the Nx curve) and, with `genome_size`, NG50/LG50. With sequences it adds the GC fraction (G+C over A+C+G+T; N and other symbols excluded), the number of N bases, the number of gaps (maximal runs of N) and other symbols. It then reports right-minus-left differences, the N50 ratio and, when both assemblies have sequences, the contigs with identical sequence on either strand (compared as sets of canonical sequences) and the shared bases. No alignment is performed, so similar but non-identical contigs are not matched.

- **API:** `dnakit.assembly.compare_assemblies(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input` (required): `[assembly_a, assembly_b]`, `{"left": a, "right": b}` (or `{"a": ..., "b": ...}`), or one assembly plus `options["other"]`; each assembly is an assembly given as FASTA text or a path to a FASTA file, `{name: sequence}`, or a list of sequences / `{"name", "sequence"}` records (contig lengths, as `{name: length}` or a list of integers, are also accepted). `options` (optional): `other` (the second assembly), `genome_size` (positive int; adds NG50 / LG50).

- **Output:** Standard envelope: `status` (`"ok"`), `module`, `operation`, `provenance`, `input_summary`, `warnings` and `result`. `result` contains `left` and `right` (statistics: `contig_count`, `total_length`, `min_length`, `max_length`, `mean_length`, `n50`, `l50`, `n90`, `l90`, `aun`, optional `ng50` / `lg50`, and with sequences `gc_fraction`, `n_bases`, `gap_count`, `other_symbol_bases`), `difference_right_minus_left`, `n50_ratio_right_to_left` and, with sequences on both sides, `shared_contig_sequences`, `shared_bases`, `left_only_contigs`, `right_only_contigs`.

- **Example code:**

```python
import json

from dnakit.assembly import compare_assemblies

left = {"c1": "AAAA", "c2": "CCCCCC"}
right = {"x": "TTTT", "y": "GGGGGGG"}
result = compare_assemblies([left, right])
print(json.dumps(result, indent=2))
```

- **Example output:**

```json
{
  "module": "assembly",
  "operation": "compare_assemblies",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 2,
    "text_length": 63,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "7b5c153e3e35be332a4aa54b930ef40bbb2ada115430172b19f7f47de8429f13"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "left": {
      "contig_count": 2,
      "total_length": 10,
      "min_length": 4,
      "max_length": 6,
      "mean_length": 5.0,
      "n50": 6,
      "l50": 1,
      "n90": 4,
      "l90": 2,
      "aun": 5.2,
      "gc_fraction": 0.6,
      "n_bases": 0,
      "gap_count": 0,
      "other_symbol_bases": 0
    },
    "right": {
      "contig_count": 2,
      "total_length": 11,
      "min_length": 4,
      "max_length": 7,
      "mean_length": 5.5,
      "n50": 7,
      "l50": 1,
      "n90": 4,
      "l90": 2,
      "aun": 5.909090909090909,
      "gc_fraction": 0.6363636363636364,
      "n_bases": 0,
      "gap_count": 0,
      "other_symbol_bases": 0
    },
    "difference_right_minus_left": {
      "contig_count": 0,
      "total_length": 1,
      "n50": 1,
      "l50": 0,
      "max_length": 1
    },
    "n50_ratio_right_to_left": 1.1666666666666667,
    "shared_contig_sequences": 1,
    "shared_bases": 4,
    "left_only_contigs": 1,
    "right_only_contigs": 1
  }
}
```

## 17) compare_assembly_graphs

- **Purpose:** Compares two GFA1 assembly graphs by their segment and link content.

- **Method:** Reads the `S` (segment) and `L` (link) lines of both GFA1 graphs. When every segment has a sequence, segments are matched by canonical sequence (the smaller of the sequence and its reverse complement); otherwise by name. Links are matched through the matched segments, a link and its reverse-complement counterpart counting as the same link; link overlaps (CIGAR) are ignored. It reports the counts, shared and unique segments and links, the Jaccard indices and whether the two graphs are identical. GFA2 and path/walk lines are not interpreted.

- **API:** `dnakit.assembly.compare_assembly_graphs(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input` (required): `[gfa_a, gfa_b]` (GFA1 text or paths) or one GFA plus `options["other"]`. `options` (optional): `other` (the second graph).

- **Output:** Standard envelope: `status` (`"ok"`), `module`, `operation`, `provenance`, `input_summary`, `warnings` and `result`. `result` contains `compared_by` (`"canonical segment sequence"` or `"segment name"`), `left` and `right` (`segments`, `links`, `total_segment_length`), `shared_segments`, `segment_jaccard`, `left_only_segments`, `right_only_segments`, `shared_links`, `link_jaccard` and `identical`.

- **Example code:**

```python
import json

from dnakit.assembly import compare_assembly_graphs

gfa_a = "S\ta\tACGT\nS\tb\tAAAC\nL\ta\t+\tb\t+\t0M\n"
gfa_b = "S\tp\tACGT\nS\tq\tGTTT\nS\tr\tCCCC\nL\tp\t+\tq\t-\t0M\nL\tq\t+\tr\t+\t0M\n"
result = compare_assembly_graphs([gfa_a, gfa_b])
print(json.dumps(result, indent=2))
```

- **Example output:**

```json
{
  "module": "assembly",
  "operation": "compare_assembly_graphs",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 2,
    "text_length": 125,
    "sequence_bases": 42,
    "lengths": [16, 26],
    "sha256": "e3c4b1e9641ea453179fce08375aa304f932df788456ec3bee4896872a1f7e64"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "compared_by": "canonical segment sequence",
    "left": {
      "segments": 2,
      "links": 1,
      "total_segment_length": 8
    },
    "right": {
      "segments": 3,
      "links": 2,
      "total_segment_length": 12
    },
    "shared_segments": 2,
    "segment_jaccard": 0.6666666666666666,
    "left_only_segments": 0,
    "right_only_segments": 1,
    "shared_links": 1,
    "link_jaccard": 0.5,
    "identical": false
  }
}
```

## 18) detect_assembly_contamination

- **Purpose:** Contamination screening is not implemented: it needs a taxonomic, vector or adapter database, or read coverage.

- **Method:** Not implemented. Deciding whether a contig is foreign requires comparison against reference databases (for example BLAST nt, Kraken2 or FCS-GX); a GC or length outlier is not a contamination call. If the contaminant contigs are already known, `remove_contamination` removes them by name. No heuristic stand-in is applied, so the function never fabricates a result; it validates nothing about the input beyond summarising it in `input_summary` and returns a structured `not_implemented` envelope.

- **API:** `dnakit.assembly.detect_assembly_contamination(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input` (required): anything (it is only summarised in `input_summary`). `options` (optional): ignored.

- **Output:** Always returns `status: "not_implemented"` (with `module`, `operation`, `provenance`, `input_summary`, `warnings`), `required_dependency` (what is missing: a taxonomic or contaminant reference database (for example BLAST nt, Kraken2, FCS-GX)) and `message` (why it cannot be computed from the input). No `result` is produced.

- **Example code:**

```python
import json

from dnakit.assembly import detect_assembly_contamination

result = detect_assembly_contamination(">c1\nACGTACGT\n>c2\nGGCCGGCC\n")
print(json.dumps(result, indent=2))
```

- **Example output:**

```json
{
  "module": "assembly",
  "operation": "detect_assembly_contamination",
  "provenance": {
    "implementation": "not_implemented",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 2,
    "text_length": 26,
    "sequence_bases": 16,
    "lengths": [8, 8],
    "sha256": "eaa8633c7014be534dfd6571d98374f8f6d24defb95b394e7f9ad5bbe9ac56f4"
  },
  "warnings": [],
  "status": "not_implemented",
  "required_dependency": "a taxonomic or contaminant reference database (for example BLAST nt, Kraken2, FCS-GX)",
  "message": "Deciding whether a contig is foreign requires comparison against reference databases; a GC or length outlier is not a contamination call."
}
```

## 19) detect_chimeric_contigs

- **Purpose:** Chimeric contig detection is not implemented: it needs read or reference alignments.

- **Method:** Not implemented. Chimeras are detected from split or discordant read alignments, or from alignments to a reference assembly, which are not part of the input. No heuristic stand-in is applied, so the function never fabricates a result; it validates nothing about the input beyond summarising it in `input_summary` and returns a structured `not_implemented` envelope.

- **API:** `dnakit.assembly.detect_chimeric_contigs(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input` (required): anything (it is only summarised in `input_summary`). `options` (optional): ignored.

- **Output:** Always returns `status: "not_implemented"` (with `module`, `operation`, `provenance`, `input_summary`, `warnings`), `required_dependency` (what is missing: read alignments (BAM/PAF) to the contigs or a reference assembly) and `message` (why it cannot be computed from the input). No `result` is produced.

- **Example code:**

```python
import json

from dnakit.assembly import detect_chimeric_contigs

result = detect_chimeric_contigs(">c1\nACGTACGT\n>c2\nGGCCGGCC\n")
print(json.dumps(result, indent=2))
```

- **Example output:**

```json
{
  "module": "assembly",
  "operation": "detect_chimeric_contigs",
  "provenance": {
    "implementation": "not_implemented",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 2,
    "text_length": 26,
    "sequence_bases": 16,
    "lengths": [8, 8],
    "sha256": "eaa8633c7014be534dfd6571d98374f8f6d24defb95b394e7f9ad5bbe9ac56f4"
  },
  "warnings": [],
  "status": "not_implemented",
  "required_dependency": "read alignments (BAM/PAF) to the contigs or a reference assembly",
  "message": "Chimeras are detected from split or discordant alignments, which are not part of the input."
}
```

## 20) detect_misassemblies

- **Purpose:** Misassembly detection is not implemented: it needs a reference alignment or read-mapping evidence.

- **Method:** Not implemented. Misassemblies are defined relative to a reference genome or to read evidence, neither of which is part of the input. No heuristic stand-in is applied, so the function never fabricates a result; it validates nothing about the input beyond summarising it in `input_summary` and returns a structured `not_implemented` envelope.

- **API:** `dnakit.assembly.detect_misassemblies(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input` (required): anything (it is only summarised in `input_summary`). `options` (optional): ignored.

- **Output:** Always returns `status: "not_implemented"` (with `module`, `operation`, `provenance`, `input_summary`, `warnings`), `required_dependency` (what is missing: a reference genome alignment (QUAST/MUMmer) or read mappings (CRAQ/Inspector)) and `message` (why it cannot be computed from the input). No `result` is produced.

- **Example code:**

```python
import json

from dnakit.assembly import detect_misassemblies

result = detect_misassemblies(">c1\nACGTACGT\n>c2\nGGCCGGCC\n")
print(json.dumps(result, indent=2))
```

- **Example output:**

```json
{
  "module": "assembly",
  "operation": "detect_misassemblies",
  "provenance": {
    "implementation": "not_implemented",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 2,
    "text_length": 26,
    "sequence_bases": 16,
    "lengths": [8, 8],
    "sha256": "eaa8633c7014be534dfd6571d98374f8f6d24defb95b394e7f9ad5bbe9ac56f4"
  },
  "warnings": [],
  "status": "not_implemented",
  "required_dependency": "a reference genome alignment (QUAST/MUMmer) or read mappings (CRAQ/Inspector)",
  "message": "Misassemblies are defined relative to a reference or to read evidence, neither of which is part of the input."
}
```

## 21) estimate_heterozygosity

- **Purpose:** Estimates the SNP heterozygosity rate of a diploid from its k-mer spectrum with a simple two-peak model.

- **Method:** Heuristic. With `peak` the homozygous coverage peak, the k-mers with coverage between `min_coverage` and `het_upper * peak` are counted as heterozygous (H, they sit near half coverage) and those in (`het_upper * peak`, `hom_upper * peak`] as homozygous (O). For isolated SNPs spaced more than k apart, H = 2*k*r*L and O = L*(1 - k*r), so with R = H/O the rate is `r = R / (k * (2 + R))`. The model ignores repeats, sequencing errors above the cutoff and clustered SNPs, and the result is only as good as the cutoff and peak; it is returned with a warning saying so. It raises an error if no k-mers fall in the homozygous window.

- **API:** `dnakit.assembly.estimate_heterozygosity(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input` (required): a k-mer histogram (`{coverage: number of distinct k-mers}`, `[[coverage, count], ...]`, or `coverage count` text/path as written by `jellyfish histo`) or sequences / reads (FASTA text, `{name: sequence}` or a list of sequences; their canonical k-mers are counted, limited to 5,000,000 bases). `options` (optional): `k` (int 3-127, default 21; needed to interpret a histogram, and used for counting sequences), `min_coverage` (error cutoff; default: the first coverage x whose count is not larger than at x+1, the first local valley) and `peak_coverage` (default: the highest bin above the cutoff, or at it when `min_coverage` is given; the lowest coverage wins ties), `het_upper` (default 0.75), `hom_upper` (default 1.5).

- **Output:** Standard envelope: `status` (`"ok"`), `module`, `operation`, `provenance`, `input_summary`, `warnings` and `result`. `result` contains `heterozygosity` (per base), `het_distinct_kmers`, `hom_distinct_kmers`, `het_to_hom_ratio`, `k`, `peak_coverage`, `min_coverage` and a `model` string; the envelope's `warnings` contains "heuristic two-peak model; see docstring for assumptions".

- **Example code:**

```python
import json

from dnakit.assembly import estimate_heterozygosity

histogram = {1: 500, 10: 20000, 20: 90000}
result = estimate_heterozygosity(histogram, options={"k": 10})
print(json.dumps(result, indent=2))
```

- **Example output:**

```json
{
  "module": "assembly",
  "operation": "estimate_heterozygosity",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 1,
    "text_length": 36,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "4c5a162b6bc453f6cdb977dc5af239f06552601c4aadcdd2717baea800cbe65b"
  },
  "warnings": ["heuristic two-peak model; see docstring for assumptions"],
  "status": "ok",
  "result": {
    "heterozygosity": 0.01,
    "het_distinct_kmers": 20000,
    "hom_distinct_kmers": 90000,
    "het_to_hom_ratio": 0.2222222222222222,
    "k": 10,
    "peak_coverage": 20,
    "min_coverage": 2,
    "model": "diploid, isolated SNPs spaced > k apart: het k-mers = 2*k*r*L, hom k-mers = L*(1-k*r); r = R / (k*(2+R)) with R = het/hom"
  }
}
```

## 22) estimate_repeats

- **Purpose:** Estimates repeat content from k-mer multiplicity: k-mers seen at more than `repeat_factor` times the unique-coverage peak count as repeats.

- **Method:** Heuristic that assumes uniform coverage. With T = sum over x >= `min_coverage` of x*n(x) (n(x) = number of distinct k-mers seen x times) and R = sum over x > `repeat_factor` * peak of x*n(x): `repeat_fraction = R/T` (share of genome positions in multi-copy k-mers), `repeat_positions_estimate = R/peak` and `genome_size_estimate = T/peak`. The estimated copy number of a repeat k-mer is `max(2, round(x/peak))`, summarised in `copy_number_classes`. Coverage variation, sequencing bias and errors above the cutoff are not modelled; the result carries a warning saying so.

- **API:** `dnakit.assembly.estimate_repeats(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input` (required): a k-mer histogram (`{coverage: number of distinct k-mers}`, `[[coverage, count], ...]`, or `coverage count` text/path as written by `jellyfish histo`) or sequences / reads (FASTA text, `{name: sequence}` or a list of sequences; their canonical k-mers are counted, limited to 5,000,000 bases). `options` (optional): `min_coverage` (error cutoff; default: the first coverage x whose count is not larger than at x+1, the first local valley) and `peak_coverage` (default: the highest bin above the cutoff, or at it when `min_coverage` is given; the lowest coverage wins ties), `repeat_factor` (> 1, default 1.5), `k` (int 3-63, default 21; only for sequence input).

- **Output:** Standard envelope: `status` (`"ok"`), `module`, `operation`, `provenance`, `input_summary`, `warnings` and `result`. `result` contains `repeat_fraction`, `repeat_positions_estimate`, `genome_size_estimate`, `repeat_distinct_kmers`, `copy_number_classes` (estimated copy number -> number of distinct k-mers), `peak_coverage`, `min_coverage`, `repeat_factor` and a `definition` string; `warnings` contains "k-mer multiplicity estimate; assumes uniform coverage".

- **Example code:**

```python
import json

from dnakit.assembly import estimate_repeats

histogram = {1: 100, 10: 1000, 20: 50, 30: 10}
result = estimate_repeats(histogram)
print(json.dumps(result, indent=2))
```

- **Example output:**

```json
{
  "module": "assembly",
  "operation": "estimate_repeats",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 1,
    "text_length": 42,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "d923a1562f061521c7247589eee78b5f69b2a419717251a52cc75a7348a63cfc"
  },
  "warnings": ["k-mer multiplicity estimate; assumes uniform coverage"],
  "status": "ok",
  "result": {
    "repeat_fraction": 0.11504424778761062,
    "repeat_positions_estimate": 130.0,
    "genome_size_estimate": 1130.0,
    "repeat_distinct_kmers": 60,
    "copy_number_classes": {
      "2": 50,
      "3": 10
    },
    "peak_coverage": 10,
    "min_coverage": 2,
    "repeat_factor": 1.5,
    "definition": "repeat_fraction = sum_{x > repeat_factor*peak} x*n(x) / sum_{x >= min_coverage} x*n(x)"
  }
}
```

## 23) estimate_size

- **Purpose:** Estimates the genome size from a k-mer histogram as the total number of k-mer instances divided by the main coverage peak.

- **Method:** Model of Liu et al. 2013: `G = sum over x >= min_coverage of x*n(x) / peak_coverage`, where n(x) is the number of distinct k-mers seen x times (optionally only up to `max_coverage`). K-mers below the cutoff are treated as sequencing errors. The result is in bases and is only as good as the cutoff and peak, which default to the first local valley and the highest bin above it.

- **API:** `dnakit.assembly.estimate_size(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input` (required): a k-mer histogram (`{coverage: number of distinct k-mers}`, `[[coverage, count], ...]`, or `coverage count` text/path as written by `jellyfish histo`) or sequences / reads (FASTA text, `{name: sequence}` or a list of sequences; their canonical k-mers are counted, limited to 5,000,000 bases). `options` (optional): `min_coverage` (error cutoff; default: the first coverage x whose count is not larger than at x+1, the first local valley) and `peak_coverage` (default: the highest bin above the cutoff, or at it when `min_coverage` is given; the lowest coverage wins ties), `max_coverage` (optional upper cutoff), `k` (int 3-63, default 21; only for sequence input).

- **Output:** Standard envelope: `status` (`"ok"`), `module`, `operation`, `provenance`, `input_summary`, `warnings` and `result`. `result` contains `genome_size` (bases), `peak_coverage`, `min_coverage`, `max_coverage`, `total_kmer_instances`, `distinct_kmers_used` and a `model` string.

- **Example code:**

```python
import json

from dnakit.assembly import estimate_size

histogram = {1: 500, 10: 20000, 20: 90000}
result = estimate_size(histogram)
print(json.dumps(result, indent=2))
```

- **Example output:**

```json
{
  "module": "assembly",
  "operation": "estimate_size",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 1,
    "text_length": 36,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "4c5a162b6bc453f6cdb977dc5af239f06552601c4aadcdd2717baea800cbe65b"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "genome_size": 100000.0,
    "peak_coverage": 20,
    "min_coverage": 2,
    "max_coverage": null,
    "total_kmer_instances": 2000000,
    "distinct_kmers_used": 110000,
    "model": "G = sum_{x>=min_coverage} x * n(x) / peak_coverage"
  }
}
```

## 24) extract_accessory_pangenome

- **Purpose:** Lists the accessory pangenome genes, i.e. those present in fewer than `core_threshold` of the genomes, split into shell and cloud genes.

- **Method:** A gene is core if it is present in at least `ceil(core_threshold * genomes)` genomes; every other gene is accessory. Accessory genes present in a fraction of genomes below `cloud_threshold` are cloud genes, the rest shell genes. Gene names are taken as already-clustered orthology groups.

- **API:** `dnakit.assembly.extract_accessory_pangenome(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input` (required): gene presence/absence data: `{genome: [genes]}`, `{gene: {genome: 0/1}}` or `{"genomes": [...], "genes": [...], "matrix": [[...]]}` (rows = genes, columns = genomes; values are booleans or non-negative numbers, > 0 meaning present). Gene names are taken as already-clustered orthology groups; no sequence clustering or homology search is done. `options` (optional): `core_threshold` (fraction in (0, 1], default 0.99), `cloud_threshold` (fraction in [0, 1], default 0.15).

- **Output:** Standard envelope: `status` (`"ok"`), `module`, `operation`, `provenance`, `input_summary`, `warnings` and `result`. `result` contains `genome_count`, `total_genes`, `core_threshold`, `accessory_genes`, `accessory_count`, `shell_genes`, `cloud_genes`, `cloud_threshold` and `gene_fractions` (fraction of genomes carrying each accessory gene).

- **Example code:**

```python
import json

from dnakit.assembly import extract_accessory_pangenome

genomes = {"G1": ["a", "b", "c"], "G2": ["a", "b", "d"], "G3": ["a", "e"]}
result = extract_accessory_pangenome(genomes, options={"core_threshold": 0.66, "cloud_threshold": 0.5})
print(json.dumps(result, indent=2))
```

- **Example output:**

```json
{
  "module": "assembly",
  "operation": "extract_accessory_pangenome",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 1,
    "text_length": 64,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "e4372679c8d0e6836c36cd361cef3caf9656e4da871974b5c7cf4f7878ed448b"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "genome_count": 3,
    "total_genes": 5,
    "core_threshold": 0.66,
    "accessory_genes": ["c", "d", "e"],
    "accessory_count": 3,
    "shell_genes": [],
    "cloud_genes": ["c", "d", "e"],
    "cloud_threshold": 0.5,
    "gene_fractions": {
      "c": 0.3333333333333333,
      "d": 0.3333333333333333,
      "e": 0.3333333333333333
    }
  }
}
```

## 25) extract_core_pangenome

- **Purpose:** Lists the core pangenome genes, i.e. those present in at least `core_threshold` of the genomes.

- **Method:** A gene is core if it is present in at least `ceil(core_threshold * genome_count)` genomes (a small tolerance of 1e-9 is applied before rounding up). With the default 0.99 a gene must be in practically all genomes; use 1.0 for strict core genes. Gene names are taken as already-clustered orthology groups.

- **API:** `dnakit.assembly.extract_core_pangenome(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input` (required): gene presence/absence data: `{genome: [genes]}`, `{gene: {genome: 0/1}}` or `{"genomes": [...], "genes": [...], "matrix": [[...]]}` (rows = genes, columns = genomes; values are booleans or non-negative numbers, > 0 meaning present). Gene names are taken as already-clustered orthology groups; no sequence clustering or homology search is done. `options` (optional): `core_threshold` (fraction in (0, 1], default 0.99).

- **Output:** Standard envelope: `status` (`"ok"`), `module`, `operation`, `provenance`, `input_summary`, `warnings` and `result`. `result` contains `genome_count`, `total_genes`, `core_threshold`, `core_min_genomes` (the minimum number of genomes), `core_genes`, `core_count` and `gene_fractions` (fraction of genomes carrying each core gene).

- **Example code:**

```python
import json

from dnakit.assembly import extract_core_pangenome

genomes = {"G1": ["a", "b", "c"], "G2": ["a", "b", "d"], "G3": ["a", "e"]}
result = extract_core_pangenome(genomes, options={"core_threshold": 0.66})
print(json.dumps(result, indent=2))
```

- **Example output:**

```json
{
  "module": "assembly",
  "operation": "extract_core_pangenome",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 1,
    "text_length": 64,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "e4372679c8d0e6836c36cd361cef3caf9656e4da871974b5c7cf4f7878ed448b"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "genome_count": 3,
    "total_genes": 5,
    "core_threshold": 0.66,
    "core_min_genomes": 2,
    "core_genes": ["a", "b"],
    "core_count": 2,
    "gene_fractions": {
      "a": 1.0,
      "b": 0.6666666666666666
    }
  }
}
```

## 26) identify_centromeres

- **Purpose:** Centromere identification is not implemented: there is no sequence-only rule, it needs satellite / repeat annotation.

- **Method:** Not implemented. Centromere positions are species-specific (satellite monomers, CENH3 ChIP, Hi-C), and no universal documented rule exists to apply to bare contigs; it needs tandem-repeat / satellite annotation (TRF, SRF, RepeatMasker) and species knowledge. No heuristic stand-in is applied, so the function never fabricates a result; it validates nothing about the input beyond summarising it in `input_summary` and returns a structured `not_implemented` envelope.

- **API:** `dnakit.assembly.identify_centromeres(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input` (required): anything (it is only summarised in `input_summary`). `options` (optional): ignored.

- **Output:** Always returns `status: "not_implemented"` (with `module`, `operation`, `provenance`, `input_summary`, `warnings`), `required_dependency` (what is missing: tandem-repeat/satellite annotation (TRF, SRF, RepeatMasker) and species-specific centromere repeat knowledge) and `message` (why it cannot be computed from the input). No `result` is produced.

- **Example code:**

```python
import json

from dnakit.assembly import identify_centromeres

result = identify_centromeres(">c1\nACGTACGT\n>c2\nGGCCGGCC\n")
print(json.dumps(result, indent=2))
```

- **Example output:**

```json
{
  "module": "assembly",
  "operation": "identify_centromeres",
  "provenance": {
    "implementation": "not_implemented",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 2,
    "text_length": 26,
    "sequence_bases": 16,
    "lengths": [8, 8],
    "sha256": "eaa8633c7014be534dfd6571d98374f8f6d24defb95b394e7f9ad5bbe9ac56f4"
  },
  "warnings": [],
  "status": "not_implemented",
  "required_dependency": "tandem-repeat/satellite annotation (TRF, SRF, RepeatMasker) and species-specific centromere repeat knowledge",
  "message": "Centromere positions are species-specific (satellite monomers, CENH3 ChIP, Hi-C); no universal documented rule exists to apply to bare contigs."
}
```

## 27) identify_repeats_in_assembly

- **Purpose:** Finds repeated regions of an assembly as the intervals covered by k-mers that occur at least `min_copies` times.

- **Method:** Counts canonical k-mers (a k-mer and its reverse complement are one) over all contigs together; k-mers containing a non-ACGT base are ignored. Every position covered by a k-mer with count >= `min_copies` is a repeat; overlapping or adjacent k-mers are merged into 0-based half-open intervals and intervals shorter than `min_length` are dropped. No repeat library and no repeat classification (transposon, satellite, ...) is used, so this is a purely internal-duplication measure; copies inside one contig and across contigs both count. Limited to 5,000,000 bases.

- **API:** `dnakit.assembly.identify_repeats_in_assembly(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input` (required): an assembly with sequences (an assembly given as FASTA text or a path to a FASTA file, `{name: sequence}`, or a list of sequences / `{"name", "sequence"}` records; lengths alone are not enough). `options` (optional): `k` (int 3-63, default 21), `min_copies` (int >= 2, default 2), `min_length` (int, default = `k`).

- **Output:** Standard envelope: `status` (`"ok"`), `module`, `operation`, `provenance`, `input_summary`, `warnings` and `result`. `result` contains `k`, `min_copies`, `min_length`, `repeat_bases`, `total_bases`, `repeat_fraction`, `repeats` (contig -> list of `{start, end, length}`) and a `definition` string.

- **Example code:**

```python
import json

from dnakit.assembly import identify_repeats_in_assembly

result = identify_repeats_in_assembly({"c": "GATTACACCGATTACA"}, options={"k": 5})
print(json.dumps(result, indent=2))
```

- **Example output:**

```json
{
  "module": "assembly",
  "operation": "identify_repeats_in_assembly",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 1,
    "text_length": 25,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "f3afd198c51e3955cb13540d57c5586de1dc62cd41af6c7fb5cc2900f62127a8"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "k": 5,
    "min_copies": 2,
    "min_length": 5,
    "repeat_bases": 14,
    "total_bases": 16,
    "repeat_fraction": 0.875,
    "repeats": {
      "c": [
        {
          "start": 0,
          "end": 7,
          "length": 7
        },
        {
          "start": 9,
          "end": 16,
          "length": 7
        }
      ]
    },
    "definition": "positions covered by canonical k-mers that occur >= min_copies times in the assembly; overlapping/adjacent k-mers merged; 0-based half-open intervals"
  }
}
```

## 28) identify_telomeres

- **Purpose:** Finds telomeric repeat arrays at contig ends as perfect tandem repeats of a telomere motif.

- **Method:** A telomere end is a perfect tandem run of at least `min_repeats` copies of the motif or of its reverse complement (any rotation of either, so `TTAGGG` also matches `CCCTAA` and `TAGGGT`) that starts within `max_distance` bases of the contig start (left end) or ends within `max_distance` bases of the contig end (right end); a run is assigned to the nearer end (ties: left), and the longest qualifying run per end is reported with 0-based half-open coordinates. Only perfect repeats are found: degenerate or short arrays are missed, and no check of chromosome-internal telomere-like repeats is made.

- **API:** `dnakit.assembly.identify_telomeres(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input` (required): an assembly with sequences (an assembly given as FASTA text or a path to a FASTA file, `{name: sequence}`, or a list of sequences / `{"name", "sequence"}` records). `options` (optional): `motif` (2-12 bases of ACGT, default `"TTAGGG"`, the vertebrate repeat; plants use `"TTTAGGG"`), `min_repeats` (int, default 5), `max_distance` (int, default 100).

- **Output:** Standard envelope: `status` (`"ok"`), `module`, `operation`, `provenance`, `input_summary`, `warnings` and `result`. `result` contains `motif`, `min_repeats`, `max_distance`, `contigs` (per contig `length`, `left` and `right` as `{start, end, copies}` or null, and `telomeres_at_ends`), `contigs_with_telomere`, `contigs_with_both_ends` and a `definition` string.

- **Example code:**

```python
import json

from dnakit.assembly import identify_telomeres

contigs = {
    "left_tel": "TTAGGG" * 6 + "ACGTACGTACGTACGGATC",
    "right_tel": "GATCGATCGATCGATC" + "CCCTAA" * 5,
    "none": "ACGTACGTACGTACGTACGTACGT",
}
result = identify_telomeres(contigs)
print(json.dumps(result, indent=2))
```

- **Example output:**

```json
{
  "module": "assembly",
  "operation": "identify_telomeres",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 1,
    "text_length": 170,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "f308faab612a4f4db560c3649114db4f013646b30b0661f9caf32b501fa001b9"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "motif": "TTAGGG",
    "min_repeats": 5,
    "max_distance": 100,
    "contigs": {
      "left_tel": {
        "length": 55,
        "left": {
          "start": 0,
          "end": 36,
          "copies": 6
        },
        "right": null,
        "telomeres_at_ends": 1
      },
      "right_tel": {
        "length": 46,
        "left": null,
        "right": {
          "start": 16,
          "end": 46,
          "copies": 5
        },
        "telomeres_at_ends": 1
      },
      "none": {
        "length": 24,
        "left": null,
        "right": null,
        "telomeres_at_ends": 0
      }
    },
    "contigs_with_telomere": 2,
    "contigs_with_both_ends": 0,
    "definition": "a telomere end is a perfect tandem run of >= min_repeats copies of the motif or its reverse complement (any rotation) lying within max_distance bases of the nearer contig end (ties go to the left end)"
  }
}
```

## 29) parse_log

- **Purpose:** Summarises an assembler log by severity level: counts of ERROR / WARNING / INFO / DEBUG lines and the first error and warning lines.

- **Method:** Each line is assigned the level of the first level word found in it (whole word, case-insensitive): CRITICAL or FATAL (counted as CRITICAL), ERROR, WARNING or WARN (counted as WARNING), INFO, DEBUG; lines without such a word are only counted in `line_count`. The first `max_lines` ERROR / CRITICAL lines and the first `max_lines` WARNING lines are returned with their line numbers (text cut at 500 characters). No tool-specific fields (SPAdes, Flye, ...) are parsed.

- **API:** `dnakit.assembly.parse_log(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input` (required): log text or a path to a log file. `options` (optional): `max_lines` (int 1-1000, default 20; how many error and how many warning lines are returned).

- **Output:** Standard envelope: `status` (`"ok"`), `module`, `operation`, `provenance`, `input_summary`, `warnings` and `result`. `result` contains `line_count`, `level_counts` (`CRITICAL`, `ERROR`, `WARNING`, `INFO`, `DEBUG`), `has_errors` (true if any ERROR or CRITICAL line), `errors` and `warnings` (lists of `{line, text}`) and a `definition` string.

- **Example code:**

```python
import json

from dnakit.assembly import parse_log

log = "INFO started\nwarning: low coverage\nERROR failed to open file\nINFO done\n[debug] details\n"
result = parse_log(log)
print(json.dumps(result, indent=2))
```

- **Example output:**

```json
{
  "module": "assembly",
  "operation": "parse_log",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 5,
    "text_length": 87,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "8e4b46eec2308c01ef1b5cd36cb9c0a55cf35e4c545db1835b4313f5c8a4cd47"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "line_count": 5,
    "level_counts": {
      "CRITICAL": 0,
      "ERROR": 1,
      "WARNING": 1,
      "INFO": 2,
      "DEBUG": 1
    },
    "has_errors": true,
    "errors": [
      {
        "line": 3,
        "text": "ERROR failed to open file"
      }
    ],
    "warnings": [
      {
        "line": 2,
        "text": "warning: low coverage"
      }
    ],
    "definition": "a line's level is the first word ERROR/CRITICAL/FATAL/WARNING/WARN/INFO/DEBUG found in it (case-insensitive)"
  }
}
```

## 30) parse_output

- **Purpose:** Parses assembler output (FASTA contigs or GFA1) into per-record length, GC and coverage fields plus assembly statistics.

- **Method:** The format is FASTA when the first non-empty line starts with `>`, otherwise GFA1. FASTA records get their length and GC fraction (G+C over A+C+G+T); SPAdes headers (`NODE_1_length_100_cov_5.2`) and MEGAHIT headers (`flag=1 multi=2.0 len=300`) are recognised and the declared length and coverage / multiplicity are returned. GFA segments get their length (from the sequence or the `LN` tag) and depth (from the `DP`, `dp` or `KC` tag, first present). Contiguity statistics (N50 etc.) are computed from the lengths. Other assemblers' header conventions are not interpreted.

- **API:** `dnakit.assembly.parse_output(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input` (required): FASTA or GFA1 text, or a path. `options`: none.

- **Output:** Standard envelope: `status` (`"ok"`), `module`, `operation`, `provenance`, `input_summary`, `warnings` and `result`. `result` contains `format` (`"fasta"` or `"gfa"`), `header_convention` (`"spades"`, `"megahit"` or null), `records` (FASTA: `name`, `length`, `gc_fraction` and, if recognised, `declared_length` with `coverage` or `multiplicity`; GFA: `name`, `length`, `depth`), `stats` (contiguity statistics as in `compare_assemblies`, from the record lengths) and, for GFA, `link_count`.

- **Example code:**

```python
import json

from dnakit.assembly import parse_output

fasta = ">NODE_1_length_8_cov_5.5\nACGTACGT\n>NODE_2_length_4_cov_2.0\nGGCC\n"
result = parse_output(fasta)
print(json.dumps(result, indent=2))
```

- **Example output:**

```json
{
  "module": "assembly",
  "operation": "parse_output",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 2,
    "text_length": 64,
    "sequence_bases": 12,
    "lengths": [8, 4],
    "sha256": "3d35bb543426e9cf6848f211f438b43042146cba719a22a02ae49eb460d310b5"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "format": "fasta",
    "header_convention": "spades",
    "records": [
      {
        "name": "NODE_1_length_8_cov_5.5",
        "length": 8,
        "gc_fraction": 0.5,
        "declared_length": 8,
        "coverage": 5.5
      },
      {
        "name": "NODE_2_length_4_cov_2.0",
        "length": 4,
        "gc_fraction": 1.0,
        "declared_length": 4,
        "coverage": 2.0
      }
    ],
    "stats": {
      "contig_count": 2,
      "total_length": 12,
      "min_length": 4,
      "max_length": 8,
      "mean_length": 6.0,
      "n50": 8,
      "l50": 1,
      "n90": 4,
      "l90": 2,
      "aun": 6.666666666666667
    }
  }
}
```

## 31) phase_assembly

- **Purpose:** Phasing an assembly is not implemented: it needs reads, Hi-C or parental data.

- **Method:** Not implemented. Assigning contigs or bases to haplotypes needs linking evidence (long reads, Hi-C or trio k-mers) and a phasing assembler such as hifiasm, HapCUT2 or WhatsHap, none of which is part of the input. No heuristic stand-in is applied, so the function never fabricates a result; it validates nothing about the input beyond summarising it in `input_summary` and returns a structured `not_implemented` envelope.

- **API:** `dnakit.assembly.phase_assembly(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input` (required): anything (it is only summarised in `input_summary`). `options` (optional): ignored.

- **Output:** Always returns `status: "not_implemented"` (with `module`, `operation`, `provenance`, `input_summary`, `warnings`), `required_dependency` (what is missing: long reads, Hi-C or trio k-mers and a phasing assembler (hifiasm, HapCUT2, WhatsHap)) and `message` (why it cannot be computed from the input). No `result` is produced.

- **Example code:**

```python
import json

from dnakit.assembly import phase_assembly

result = phase_assembly(">c1\nACGTACGT\n>c2\nGGCCGGCC\n")
print(json.dumps(result, indent=2))
```

- **Example output:**

```json
{
  "module": "assembly",
  "operation": "phase_assembly",
  "provenance": {
    "implementation": "not_implemented",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 2,
    "text_length": 26,
    "sequence_bases": 16,
    "lengths": [8, 8],
    "sha256": "eaa8633c7014be534dfd6571d98374f8f6d24defb95b394e7f9ad5bbe9ac56f4"
  },
  "warnings": [],
  "status": "not_implemented",
  "required_dependency": "long reads, Hi-C or trio k-mers and a phasing assembler (hifiasm, HapCUT2, WhatsHap)",
  "message": "Assigning contigs or bases to haplotypes needs linking evidence that is not part of the input."
}
```

## 32) plot_graph

- **Purpose:** Draws a GFA1 assembly graph as an SVG with the segments on a circle and the links as straight lines.

- **Method:** Segments are placed at equal angles on a circle in sorted name order (640 x 640 px); node radius is `4 + 2*log10(length + 1)` pixels; each link is a straight line between its two segments (orientation is not drawn); names are written next to the nodes only when there are at most 60 segments, and every node has a tooltip with its name and length. The layout is deterministic and carries no biological meaning (it is not a Bandage-style layout). At most 2,000 segments; pure Python. If `output` is given the SVG is also written to that file.

- **API:** `dnakit.assembly.plot_graph(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input` (required): GFA1 text or a path (at most 2,000 segments). `options` (optional): `output` (file path to write the SVG; default none).

- **Output:** Standard envelope: `status` (`"ok"`), `module`, `operation`, `provenance`, `input_summary`, `warnings` and `result`. `result` contains `svg` (the SVG text), `node_count`, `edge_count` and `output_path` (the file written, or null).

- **Example code:**

```python
import json

from dnakit.assembly import plot_graph

gfa = "S\tp\tACGT\nS\tq\tGTTT\nS\tr\tCCCC\nL\tp\t+\tq\t-\t0M\nL\tq\t+\tr\t+\t0M\n"
result = plot_graph(gfa)
print(json.dumps(result, indent=2))
```

- **Example output:**

```json
{
  "module": "assembly",
  "operation": "plot_graph",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 5,
    "text_length": 53,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "32ffc4c7e7ccc5721dfef1ca795811be545710992f65d6e40a4a80807358f2e0"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "svg": "<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"640\" height=\"640\" viewBox=\"0 0 640 640\">\n<line x1=\"600.0\" y1=\"320.0\" x2=\"180.0\" y2=\"562.5\" stroke=\"#888\" stroke-width=\"1\"/>\n<line x1=\"180.0\" y1=\"562.5\" x2=\"180.0\" y2=\"77.5\" stroke=\"#888\" stroke...",
    "node_count": 3,
    "edge_count": 2,
    "output_path": null
  }
}
```

(Long arrays and text in the example output are shortened with `...`.)

## 33) polish

- **Purpose:** Polishing is not implemented: it needs reads aligned to the assembly and a polisher such as Racon, Medaka or Pilon.

- **Method:** Not implemented. Consensus correction is computed from read pileups over the assembly, which are not part of the input. No heuristic stand-in is applied, so the function never fabricates a result; it validates nothing about the input beyond summarising it in `input_summary` and returns a structured `not_implemented` envelope.

- **API:** `dnakit.assembly.polish(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input` (required): anything (it is only summarised in `input_summary`). `options` (optional): ignored.

- **Output:** Always returns `status: "not_implemented"` (with `module`, `operation`, `provenance`, `input_summary`, `warnings`), `required_dependency` (what is missing: sequencing reads, an aligner and a polishing tool (Racon, Medaka, Pilon, NextPolish)) and `message` (why it cannot be computed from the input). No `result` is produced.

- **Example code:**

```python
import json

from dnakit.assembly import polish

result = polish(">c1\nACGTACGT\n>c2\nGGCCGGCC\n")
print(json.dumps(result, indent=2))
```

- **Example output:**

```json
{
  "module": "assembly",
  "operation": "polish",
  "provenance": {
    "implementation": "not_implemented",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 2,
    "text_length": 26,
    "sequence_bases": 16,
    "lengths": [8, 8],
    "sha256": "eaa8633c7014be534dfd6571d98374f8f6d24defb95b394e7f9ad5bbe9ac56f4"
  },
  "warnings": [],
  "status": "not_implemented",
  "required_dependency": "sequencing reads, an aligner and a polishing tool (Racon, Medaka, Pilon, NextPolish)",
  "message": "Consensus correction is computed from read pileups, which are not part of the input."
}
```

## 34) purge_haplotigs

- **Purpose:** Purging haplotigs is not implemented: it needs read-depth and self-alignment evidence (purge_dups, Purge_Haplotigs).

- **Method:** Not implemented. Haplotigs are identified from half-coverage and allelic self-alignment evidence that is not part of the input. No heuristic stand-in is applied, so the function never fabricates a result; it validates nothing about the input beyond summarising it in `input_summary` and returns a structured `not_implemented` envelope.

- **API:** `dnakit.assembly.purge_haplotigs(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input` (required): anything (it is only summarised in `input_summary`). `options` (optional): ignored.

- **Output:** Always returns `status: "not_implemented"` (with `module`, `operation`, `provenance`, `input_summary`, `warnings`), `required_dependency` (what is missing: read coverage and contig self-alignments (purge_dups, Purge_Haplotigs)) and `message` (why it cannot be computed from the input). No `result` is produced.

- **Example code:**

```python
import json

from dnakit.assembly import purge_haplotigs

result = purge_haplotigs(">c1\nACGTACGT\n>c2\nGGCCGGCC\n")
print(json.dumps(result, indent=2))
```

- **Example output:**

```json
{
  "module": "assembly",
  "operation": "purge_haplotigs",
  "provenance": {
    "implementation": "not_implemented",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 2,
    "text_length": 26,
    "sequence_bases": 16,
    "lengths": [8, 8],
    "sha256": "eaa8633c7014be534dfd6571d98374f8f6d24defb95b394e7f9ad5bbe9ac56f4"
  },
  "warnings": [],
  "status": "not_implemented",
  "required_dependency": "read coverage and contig self-alignments (purge_dups, Purge_Haplotigs)",
  "message": "Haplotigs are identified from half-coverage and allelic alignment evidence that is not part of the input."
}
```

## 35) remove_contamination

- **Purpose:** Removes the contigs named in `options.contaminants` from an assembly (detection is not performed).

- **Method:** Drops the contigs whose names are in the list and keeps the rest in their original order; names that are not in the assembly are listed in `not_found` rather than raising an error. It reports the contiguity statistics before and after, and the FASTA text of the kept contigs when sequences are available. Removing every contig raises an error. The list of contaminants must come from a screen done elsewhere (`detect_assembly_contamination` is not implemented).

- **API:** `dnakit.assembly.remove_contamination(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input` (required): an assembly given as FASTA text or a path to a FASTA file, `{name: sequence}`, or a list of sequences / `{"name", "sequence"}` records (contig lengths, as `{name: length}` or a list of integers, are also accepted). `options` (required key): `contaminants` (list of contig names to remove).

- **Output:** Standard envelope: `status` (`"ok"`), `module`, `operation`, `provenance`, `input_summary`, `warnings` and `result`. `result` contains `kept` and `removed` (contig names), `not_found`, `stats_before`, `stats_after` (contiguity statistics as in `compare_assemblies`) and, when sequences are available, `fasta` (the kept contigs).

- **Example code:**

```python
import json

from dnakit.assembly import remove_contamination

assembly = {"c1": "AAAA", "c2": "CCCCCC", "c3": "GG"}
result = remove_contamination(assembly, options={"contaminants": ["c2", "zz"]})
print(json.dumps(result, indent=2))
```

- **Example output:**

```json
{
  "module": "assembly",
  "operation": "remove_contamination",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 1,
    "text_length": 42,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "5a945aa9f35559b3a7f890ec3868e4c370d32bbc5d5dbdce9aa4c71273e9bb14"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "kept": ["c1", "c3"],
    "removed": ["c2"],
    "not_found": ["zz"],
    "stats_before": {
      "contig_count": 3,
      "total_length": 12,
      "min_length": 2,
      "max_length": 6,
      "mean_length": 4.0,
      "n50": 6,
      "l50": 1,
      "n90": 2,
      "l90": 3,
      "aun": 4.666666666666667,
      "gc_fraction": 0.6666666666666666,
      "n_bases": 0,
      "gap_count": 0,
      "other_symbol_bases": 0
    },
    "stats_after": {
      "contig_count": 2,
      "total_length": 6,
      "min_length": 2,
      "max_length": 4,
      "mean_length": 3.0,
      "n50": 4,
      "l50": 1,
      "n90": 2,
      "l90": 2,
      "aun": 3.3333333333333335,
      "gc_fraction": 0.3333333333333333,
      "n_bases": 0,
      "gap_count": 0,
      "other_symbol_bases": 0
    },
    "fasta": ">c1\nAAAA\n>c3\nGG\n"
  }
}
```

## 36) run

- **Purpose:** Runs one caller-installed assembler (Flye, hifiasm, MEGAHIT, Shasta or Unicycler) through the matching `assemble_*` wrapper.

- **Method:** First validates the settings exactly as `build_command` does (so a missing required key raises an error even when the program is absent), then looks the executable up on PATH. If it is not found the function returns `not_implemented` naming the program and nothing is run. If it is found, the matching wrapper (`assemble_flye`, `assemble_hifiasm`, `assemble_megahit`, `assemble_shasta` or `assemble_unicycler`) is called with the same argument vector: no shell, output only under `output_directory`, with a timeout (`timeout_seconds`) and output limits; its result is returned under `result`. DNAKit does not install or download any assembler. The example below was produced without any assembler installed, so it shows the `not_implemented` path (the success path was exercised in the unit tests with a stand-in executable, not with a real assembler).

- **API:** `dnakit.assembly.run(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input` (required): a mapping of settings (`options` are merged over it and win): `tool` (`flye`, `hifiasm`, `megahit`, `shasta` or `unicycler`), `output_directory`, and the tool's inputs and parameters exactly as in `build_command` (`reads`, `read1`, `read2`, `long_reads`, `genome_size`, `platform`, `config`, `mode`, `prefix`, `threads`, `min_contig_length`, `program`), plus an optional `timeout_seconds`.

- **Output:** Without the executable on PATH: `status: "not_implemented"` (with `module`, `operation`, `provenance`, `input_summary`, `warnings`), `required_dependency` (`"<program> executable on PATH"`) and `message` (the program is not installed or not on PATH, so nothing was run); no `result`. With the executable: the standard envelope with `status: "ok"` and `result` holding the wrapper's report (`status`, `operation` such as `assemble_flye`, `backend`, a `result` with the output directory, the assembly files, `artifacts` and the program `log`, and `provenance`).

- **Example code:**

```python
import json

from dnakit.assembly import run

settings = {
    "tool": "flye",
    "reads": "reads.fastq",
    "genome_size": "5m",
    "output_directory": "flye_out",
}
result = run(settings)
print(json.dumps(result, indent=2))
```

- **Example output:**

```json
{
  "module": "assembly",
  "operation": "run",
  "provenance": {
    "implementation": "not_implemented",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 1,
    "text_length": 93,
    "sequence_bases": 10,
    "lengths": [10],
    "sha256": "858ba72be547f61d29e0e8be6d87531c90bc5c6b4be8918a4ca7c84b978670fd"
  },
  "warnings": [],
  "status": "not_implemented",
  "required_dependency": "flye executable on PATH",
  "message": "The assembler program 'flye' is not installed or not on PATH, so nothing was run."
}
```

## 37) scaffold

- **Purpose:** Joins contigs into scaffolds following an explicit, caller-supplied layout (no linking evidence is inferred).

- **Method:** For each scaffold in `layout` the listed contigs are concatenated in order; a contig with orientation `-` is reverse-complemented first, and consecutive contigs are separated by a run of `N` of the given length (per-contig `gap`, else the global `gap`; a gap of 0 inserts nothing). A contig may be used only once. Ordering and orientation are NOT predicted: the layout must come from other evidence (for example a Hi-C scaffolder). The AGP-like table uses 0-based half-open coordinates on each scaffold.

- **API:** `dnakit.assembly.scaffold(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input` (required): an assembly with sequences (an assembly given as FASTA text or a path to a FASTA file, `{name: sequence}`, or a list of sequences / `{"name", "sequence"}` records). `options`: `layout` (required): a list of scaffolds `{"name": str, "contigs": [contig_name | {"contig": name, "orientation": "+" | "-", "gap": int}, ...]}` where `gap` is the number of N between this contig and the next; `gap` (int >= 0, default 100).

- **Output:** Standard envelope: `status` (`"ok"`), `module`, `operation`, `provenance`, `input_summary`, `warnings` and `result`. `result` contains `scaffolds` (each with `name`, `sequence`, `length`, `contig_count`), `agp` (rows with `scaffold`, `start`, `end`, `type` `contig` or `gap`, `contig`, `orientation`) and `unplaced_contigs` (contigs not used in the layout).

- **Example code:**

```python
import json

from dnakit.assembly import scaffold

contigs = {"c1": "AACC", "c2": "GGTT", "c3": "ACGT"}
layout = [{"name": "S1", "contigs": ["c1", {"contig": "c2", "orientation": "-"}]}]
result = scaffold(contigs, options={"layout": layout, "gap": 5})
print(json.dumps(result, indent=2))
```

- **Example output:**

```json
{
  "module": "assembly",
  "operation": "scaffold",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 1,
    "text_length": 42,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "6646b5dc42423ea3aa5dd2d01c9954cd36aa9d53434aa8fc7dafd4747ac80aa4"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "scaffolds": [
      {
        "name": "S1",
        "sequence": "AACCNNNNNAACC",
        "length": 13,
        "contig_count": 2
      }
    ],
    "agp": [
      {
        "scaffold": "S1",
        "start": 0,
        "end": 4,
        "type": "contig",
        "contig": "c1",
        "orientation": "+"
      },
      {
        "scaffold": "S1",
        "start": 4,
        "end": 9,
        "type": "gap",
        "contig": null,
        "orientation": null
      },
      {
        "scaffold": "S1",
        "start": 9,
        "end": 13,
        "type": "contig",
        "contig": "c2",
        "orientation": "-"
      }
    ],
    "unplaced_contigs": ["c3"]
  }
}
```
