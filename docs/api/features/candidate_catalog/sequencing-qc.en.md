# Sequencing quality control

This category contains **50** documented capabilities. 46 are native DNAKit implementations in pure Python: 42 compute a result from the supplied reads, SAM text, fragment tables or counts and return `status: "ok"`, and 4 (`basecall_nanopore_reads`, `calculate_phasing_qc`, `calculate_phasing_qc_tenx`, `correct_long_reads`) always return `status: "not_implemented"` with `required_dependency` and `message` because they need a basecaller, Illumina InterOp files or a long-read correction program. The other 4 (`run_cutadapt`, `run_fastp`, `run_fastqc`, `run_multiqc`) are bounded wrappers that call a caller-installed external program. The native metrics follow FastQC, Picard, umi_tools and Cell Ranger definitions on the supplied reads only (no alignment), and several are labelled heuristics in their sections. The functions are exported through the formal DNAKit namespace and also exposed through the Agent/MCP registry where the JSON call plan is compatible. The example outputs of the native functions are real outputs; long arrays and strings are abbreviated with `...`.

## 1) assign_sample_barcodes

- **Purpose:** Assigns each read to a sample by matching its barcode against a sample sheet with a mismatch tolerance.

- **Method:** The barcode is taken from the read start, the read end or the Illumina header comment (the index at the end of the comment, `i7+i5` for dual indexes) and compared with every sample barcode by Hamming distance, with `N` counting as a mismatch. A read is assigned to the closest sample if that distance is at most `max_mismatches` and the best sample is closer than the second best by at least `min_distance_gap` mismatches; reads with no barcode within tolerance are `unassigned` and reads that fail the gap rule are `ambiguous`. Pure Python.

- **API:** `dnakit.sequencing_qc.assign_sample_barcodes(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input` is FASTQ text or a file path (`.gz` allowed), a list of reads ((id, sequence, quality) tuples, `{sequence, quality}` mappings, FASTQ record objects or bare sequences). A pair mapping is not accepted (`INVALID_INPUT`); for paired data pass the read that carries the sample barcode (usually the index read or read 1). Options: `samples` (required; mapping sample -> barcode, dual indexes written `i7+i5`), `barcode_location` (`read_start` default, `read_end` or `header`), `barcode_length` (default: the sample-sheet barcode length), `max_mismatches` (default 1), `min_distance_gap` (default 1), `max_listed` (per-read rows returned, default 100000), `phred_offset` (33 or 64, default 33), `max_reads` (default 1000000).

- **Output:** Returns a JSON-serializable envelope with `status` (`"ok"`), `module`, `operation`, `provenance` (`implementation: "native"`, `backend: "local"`), `input_summary`, `warnings` and `result`. `result` holds `reads`, `barcode_location`, `max_mismatches`, `min_distance_gap`, `reads_per_sample`, `assigned`, `unassigned`, `ambiguous`, `assigned_fraction`, `min_distance_between_sample_barcodes` and `assignments` (per read `id`, `barcode`, `sample`, `distance`, `status`). Invalid input or options raise `InputFormatError` or `ConfigurationError`.

- **Example code:**

```python
from dnakit.sequencing_qc import assign_sample_barcodes

def fastq(*records):
    return "".join(f"@{name}\n{seq}\n+\n{qual}\n" for name, seq, qual in records)

reads = fastq(
    ("a", "ACGTAAAA", "IIIIJJJJ"),
    ("b", "ACGAAAAA", "IIIIJJJJ"),
    ("c", "TTTTCCCC", "IIIIJJJJ"),
    ("d", "GGGGCCCC", "IIIIJJJJ"),
)
result = assign_sample_barcodes(reads, options={"samples": {"S1": "ACGT", "S2": "TTTT"}, "barcode_length": 4})
print(result)
```

- **Example output:**

```json
{
  "module": "sequencing_qc",
  "operation": "assign_sample_barcodes",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 4,
    "text_length": 92,
    "sequence_bases": 32,
    "lengths": [8, 8, 8, 8],
    "sha256": "c5f05df413cd46a9d8fecf917a84d8586235800f704f7497bc3ea55d1dfcc2ad"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "reads": 4,
    "barcode_location": "read_start",
    "max_mismatches": 1,
    "min_distance_gap": 1,
    "reads_per_sample": {
      "S1": 2,
      "S2": 1
    },
    "assigned": 3,
    "unassigned": 1,
    "ambiguous": 0,
    "assigned_fraction": 0.75,
    "min_distance_between_sample_barcodes": 3,
    "assignments": [
      {
        "id": "a",
        "barcode": "ACGT",
        "sample": "S1",
        "distance": 0,
        "status": "assigned"
      },
      {
        "id": "b",
        "barcode": "ACGA",
        "sample": "S1",
        "distance": 1,
        "status": "assigned"
      },
      "..."
    ]
  }
}
```

## 2) basecall_nanopore_reads

- **Purpose:** Intended to basecall Oxford Nanopore raw signal (FAST5/POD5) into bases with a neural-network basecaller.

- **Method:** Not implemented. Basecalling converts raw current signal to bases with a neural-network model, and no basecaller or model is bundled with DNAKit. The function always returns `status: "not_implemented"` and performs no computation.

- **API:** `dnakit.sequencing_qc.basecall_nanopore_reads(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input` is accepted but not used for any computation (a path or any object); `options` is optional and ignored.

- **Output:** Always returns an envelope with `status: "not_implemented"`, `module`, `operation`, `provenance` (`implementation: "not_implemented"`), `input_summary`, `warnings`, `required_dependency` ("Oxford Nanopore basecaller (dorado/guppy) and trained model weights") and `message` (the reason). There is no `result` key.

- **Example code:**

```python
from dnakit.sequencing_qc import basecall_nanopore_reads

result = basecall_nanopore_reads("run1.pod5")
print(result)
```

- **Example output:**

```json
{
  "module": "sequencing_qc",
  "operation": "basecall_nanopore_reads",
  "provenance": {
    "implementation": "not_implemented",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 1,
    "text_length": 9,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "10557ad64d8a6c03ccd022c14b535ba4280940dd255f1a2f5cede803ce64fa27"
  },
  "warnings": [],
  "status": "not_implemented",
  "required_dependency": "Oxford Nanopore basecaller (dorado/guppy) and trained model weights",
  "message": "Basecalling converts raw FAST5/POD5 current signal to bases with a neural-network model; no basecaller or model is bundled."
}
```

## 3) calculate_adapter_content

- **Purpose:** Calculates the cumulative adapter content along reads (FastQC-style k-mer rule).

- **Method:** A read hits an adapter when the first min(`kmer_size`, adapter length) bases of the adapter occur in the read (Hamming distance at most `max_mismatches`). The cumulative percentage at position p counts the reads whose first hit starts at or before p (1-based), as in FastQC. Adapters truncated by the read end are not detected. Pure Python.

- **API:** `dnakit.sequencing_qc.calculate_adapter_content(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input` is FASTQ text or a file path (`.gz` allowed), a list of reads ((id, sequence, quality) tuples, `{sequence, quality}` mappings, FASTQ record objects or bare sequences), or a pair `{"read1": ..., "read2": ...}` (analysed per mate). Options: `adapters` (mapping name -> sequence or list of sequences; default FastQC's Illumina universal, small-RNA 3' and 5', Nextera and SOLiD small-RNA adapters), `kmer_size` (default 12), `max_mismatches` (default 0), `phred_offset` (33 or 64, default 33), `max_reads` (default 1000000).

- **Output:** Returns a JSON-serializable envelope with `status` (`"ok"`), `module`, `operation`, `provenance` (`implementation: "native"`, `backend: "local"`), `input_summary`, `warnings` and `result`. `result` holds `reads`, `kmer_size`, `max_mismatches`, `reads_with_any_adapter`, `fraction_reads_with_any_adapter`, `adapters` (per adapter `sequence`, `kmer_used`, `reads_with_adapter`, `fraction_reads_with_adapter`, `max_cumulative_percent`, `cumulative_percent_by_position`) and a `definition` string. For a pair the result has `paired: true` and holds these fields separately under `read1` and `read2`. Invalid input or options raise `InputFormatError` or `ConfigurationError`.

- **Example code:**

```python
from dnakit.sequencing_qc import calculate_adapter_content

reads = ["AAAAAGATCGGAAGAGCTT", "A" * 19]
result = calculate_adapter_content(reads, options={"adapters": {"Illumina Universal Adapter": "AGATCGGAAGAG"}})
print(result)
```

- **Example output:**

```json
{
  "module": "sequencing_qc",
  "operation": "calculate_adapter_content",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 2,
    "text_length": 46,
    "sequence_bases": 38,
    "lengths": [
      19,
      19
    ],
    "sha256": "fbd53e50eb19d3d8570f3ff850b3a816ffdeb3ed41eb8f31c591d1957cee9b9c"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "reads": 2,
    "kmer_size": 12,
    "max_mismatches": 0,
    "reads_with_any_adapter": 1,
    "fraction_reads_with_any_adapter": 0.5,
    "adapters": {
      "Illumina Universal Adapter": {
        "sequence": "AGATCGGAAGAG",
        "kmer_used": "AGATCGGAAGAG",
        "reads_with_adapter": 1,
        "fraction_reads_with_adapter": 0.5,
        "max_cumulative_percent": 50.0,
        "cumulative_percent_by_position": [
          0.0,
          0.0,
          "..."
        ]
      }
    },
    "definition": "A read hits an adapter when the first min(kmer_size, len(adapter)) bases of the adapter occur in the read (Hamming <= max_mismatches). The cumulative percentage..."
  }
}
```

## 4) calculate_adapter_content_tenx

- **Purpose:** Calculates adapter, template-switch-oligo and poly-A content of 10x Chromium read 2 (cDNA).

- **Method:** Read 2 is checked with the same k-mer rule as `calculate_adapter_content`; the default adapter list adds the 10x template switch oligo (AAGCAGTGGTATCAACGCAGAGTACATGGG) and a 12-mer poly-A. The 10x chemistry is not auto-detected: `chemistry` (or the explicit lengths) decides how read 1 is split into cell barcode and UMI.

- **API:** `dnakit.sequencing_qc.calculate_adapter_content_tenx(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input` is a pair `{"read1": <barcode+UMI reads>, "read2": <cDNA reads>}` (FASTQ text, path or read list); a single read set is rejected and read 1 must be at least `barcode_length` + `umi_length` long. Options: `chemistry` (`v1` 14+10, `v2` 16+10, `v3` 16+12 bp barcode+UMI; default `v3`), `barcode_length`, `umi_length`, `phred_offset` (33), `max_reads` (1000000). Extra options as `calculate_adapter_content`: `adapters`, `kmer_size` (12), `max_mismatches` (0).

- **Output:** Returns a JSON-serializable envelope with `status` (`"ok"`), `module`, `operation`, `provenance` (`implementation: "native"`, `backend: "local"`), `input_summary`, `warnings` and `result`. `result` holds `chemistry`, `barcode_length`, `umi_length`, `read_pairs`, `region` (`rna_read (read 2)`) and the fields of `calculate_adapter_content` for read 2 (`reads`, `kmer_size`, `max_mismatches`, `reads_with_any_adapter`, `fraction_reads_with_any_adapter`, `adapters`, `definition`). Invalid input or options raise `InputFormatError` or `ConfigurationError`.

- **Example code:**

```python
from dnakit.sequencing_qc import calculate_adapter_content_tenx

def fastq(*records):
    return "".join(f"@{name}\n{seq}\n+\n{qual}\n" for name, seq, qual in records)

barcode = "ACGTACGTACGTACGT"
tenx = {
    "read1": fastq(
        ("p1", barcode + "A" * 12, "I" * 28),
        ("p2", barcode + "A" * 12, "I" * 28),
        ("p3", "T" * 16 + "C" * 12, "I" * 16 + "5" * 12),
    ),
    "read2": fastq(
        ("p1", "GATTACAGATCGGAAGAGCGGTTCAGCAGG", "I" * 30),
        ("p2", "GATTACAGATCGGAAGAGCGGTTCAGCAGG", "I" * 30),
        ("p3", "C" * 30, "5" * 30),
    ),
}
result = calculate_adapter_content_tenx(tenx, options={"adapters": {"Illumina Universal Adapter": "AGATCGGAAGAG"}})
print(result)
```

- **Example output:**

```json
{
  "module": "sequencing_qc",
  "operation": "calculate_adapter_content_tenx",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 1,
    "text_length": 446,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "c2bd8879256656753e784d4d4fe4f7d2de6add20603c60c539095f494a8d364c"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "chemistry": "v3",
    "barcode_length": 16,
    "umi_length": 12,
    "read_pairs": 3,
    "region": "rna_read (read 2)",
    "reads": 3,
    "kmer_size": 12,
    "max_mismatches": 0,
    "reads_with_any_adapter": 2,
    "fraction_reads_with_any_adapter": 0.666667,
    "adapters": {
      "Illumina Universal Adapter": {
        "sequence": "AGATCGGAAGAG",
        "kmer_used": "AGATCGGAAGAG",
        "reads_with_adapter": 2,
        "fraction_reads_with_adapter": 0.666667,
        "max_cumulative_percent": 66.666667,
        "cumulative_percent_by_position": [
          0.0,
          0.0,
          "..."
        ]
      }
    },
    "definition": "A read hits an adapter when the first min(kmer_size, len(adapter)) bases of the adapter occur in the read (Hamming <= max_mismatches). The cumulative percentage..."
  }
}
```

## 5) calculate_barcode_quality

- **Purpose:** Calculates the sequencing quality and validity of read barcodes (sample or cell barcodes).

- **Method:** The barcode is cut from the read start, the read end or the header comment. Reported are the number of distinct barcodes, those containing `N`, the most frequent barcodes, the Shannon entropy of the barcode counts in bits and, when qualities exist, the mean quality, Q30 fraction and mean quality per position. With a `whitelist` each barcode is classified as an exact match, corrected (a unique whitelist barcode within `max_mismatches`), ambiguous (several equally close) or invalid. Pure Python.

- **API:** `dnakit.sequencing_qc.calculate_barcode_quality(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input` is FASTQ text or a file path (`.gz` allowed), a list of reads ((id, sequence, quality) tuples, `{sequence, quality}` mappings, FASTQ record objects or bare sequences). A pair mapping is not accepted (`INVALID_INPUT`); for paired data pass the read that carries the barcode. Options: `barcode_location` (`read_start` default, `read_end` or `header`), `barcode_length` (default 8; not used for `header`), `whitelist` (list or newline-separated text of valid barcodes), `max_mismatches` (default 1, used with a whitelist), `phred_offset` (33 or 64, default 33), `max_reads` (default 1000000).

- **Output:** Returns a JSON-serializable envelope with `status` (`"ok"`), `module`, `operation`, `provenance` (`implementation: "native"`, `backend: "local"`), `input_summary`, `warnings` and `result`. `result` holds `reads`, `barcode_source`, `distinct_barcodes`, `reads_with_N_in_barcode`, `top_barcodes`, `barcode_entropy_bits`, `mean_quality`, `q30_fraction`, `mean_quality_by_position` and, with a whitelist, `whitelist` (`size`, `max_mismatches`, `exact_match_fraction`, `corrected_fraction`, `ambiguous_fraction`, `invalid_fraction`). Invalid input or options raise `InputFormatError` or `ConfigurationError`.

- **Example code:**

```python
from dnakit.sequencing_qc import calculate_barcode_quality

def fastq(*records):
    return "".join(f"@{name}\n{seq}\n+\n{qual}\n" for name, seq, qual in records)

reads = fastq(
    ("a", "ACGTTTTT", "IIII####"),
    ("b", "ACGATTTT", "IIII####"),
    ("c", "TTTTAAAA", "IIII####"),
    ("d", "ACGTCCCC", "IIII####"),
)
result = calculate_barcode_quality(reads, options={"barcode_length": 4, "whitelist": ["ACGT", "TTTT"]})
print(result)
```

- **Example output:**

```json
{
  "module": "sequencing_qc",
  "operation": "calculate_barcode_quality",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 4,
    "text_length": 92,
    "sequence_bases": 32,
    "lengths": [8, 8, 8, 8],
    "sha256": "1ea787f978a89222ad3f525a54d47607e60402c92e5b3e1e32f8f3633b41e325"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "reads": 4,
    "barcode_source": "read_start",
    "distinct_barcodes": 3,
    "reads_with_N_in_barcode": 0,
    "top_barcodes": [
      {
        "barcode": "ACGT",
        "reads": 2
      },
      {
        "barcode": "ACGA",
        "reads": 1
      },
      "..."
    ],
    "barcode_entropy_bits": 1.5,
    "mean_quality": 40.0,
    "q30_fraction": 1.0,
    "mean_quality_by_position": [
      40.0,
      40.0,
      "..."
    ],
    "whitelist": {
      "size": 2,
      "max_mismatches": 1,
      "exact_match_fraction": 0.75,
      "corrected_fraction": 0.25,
      "ambiguous_fraction": 0.0,
      "invalid_fraction": 0.0
    }
  }
}
```

## 6) calculate_barcode_quality_tenx

- **Purpose:** Calculates quality and whitelist validity of 10x cell barcodes and UMIs from read 1.

- **Method:** Read 1 is split into the cell barcode and the UMI by the chemistry layout. The cell barcode is reported as in `calculate_barcode_quality` (with the optional whitelist); the UMI gets its distinct count, `N` count, mean quality, Q30 fraction and the number of reads with a homopolymer UMI. The 10x chemistry is not auto-detected: `chemistry` (or the explicit lengths) decides how read 1 is split into cell barcode and UMI.

- **API:** `dnakit.sequencing_qc.calculate_barcode_quality_tenx(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input` is a pair `{"read1": <barcode+UMI reads>, "read2": <cDNA reads>}` (FASTQ text, path or read list); a single read set is rejected and read 1 must be at least `barcode_length` + `umi_length` long. Options: `chemistry` (`v1` 14+10, `v2` 16+10, `v3` 16+12 bp barcode+UMI; default `v3`), `barcode_length`, `umi_length`, `phred_offset` (33), `max_reads` (1000000). Extra options: `whitelist` (valid cell barcodes), `max_mismatches` (default 1).

- **Output:** Returns a JSON-serializable envelope with `status` (`"ok"`), `module`, `operation`, `provenance` (`implementation: "native"`, `backend: "local"`), `input_summary`, `warnings` and `result`. `result` holds `chemistry`, `barcode_length`, `umi_length`, `read_pairs`, `cell_barcode` (the report of `calculate_barcode_quality`) and `umi` (`distinct_umis`, `reads_with_N_in_umi`, `mean_quality`, `q30_fraction`, `homopolymer_umi_reads`). Invalid input or options raise `InputFormatError` or `ConfigurationError`.

- **Example code:**

```python
from dnakit.sequencing_qc import calculate_barcode_quality_tenx

def fastq(*records):
    return "".join(f"@{name}\n{seq}\n+\n{qual}\n" for name, seq, qual in records)

barcode = "ACGTACGTACGTACGT"
tenx = {
    "read1": fastq(
        ("p1", barcode + "A" * 12, "I" * 28),
        ("p2", barcode + "A" * 12, "I" * 28),
        ("p3", "T" * 16 + "C" * 12, "I" * 16 + "5" * 12),
    ),
    "read2": fastq(
        ("p1", "GATTACAGATCGGAAGAGCGGTTCAGCAGG", "I" * 30),
        ("p2", "GATTACAGATCGGAAGAGCGGTTCAGCAGG", "I" * 30),
        ("p3", "C" * 30, "5" * 30),
    ),
}
result = calculate_barcode_quality_tenx(tenx, options={"whitelist": [barcode]})
print(result)
```

- **Example output:**

```json
{
  "module": "sequencing_qc",
  "operation": "calculate_barcode_quality_tenx",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 1,
    "text_length": 446,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "c2bd8879256656753e784d4d4fe4f7d2de6add20603c60c539095f494a8d364c"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "chemistry": "v3",
    "barcode_length": 16,
    "umi_length": 12,
    "read_pairs": 3,
    "cell_barcode": {
      "reads": 3,
      "barcode_source": "read_start",
      "distinct_barcodes": 2,
      "reads_with_N_in_barcode": 0,
      "top_barcodes": [
        {
          "barcode": "ACGTACGTACGTACGT",
          "reads": 2
        },
        {
          "barcode": "TTTTTTTTTTTTTTTT",
          "reads": 1
        }
      ],
      "barcode_entropy_bits": 0.918296,
      "mean_quality": 40.0,
      "q30_fraction": 1.0,
      "mean_quality_by_position": [
        40.0,
        40.0,
        "..."
      ],
      "whitelist": {
        "size": 1,
        "max_mismatches": 1,
        "exact_match_fraction": 0.666667,
        "corrected_fraction": 0.0,
        "ambiguous_fraction": 0.0,
        "invalid_fraction": 0.333333
      }
    },
    "umi": {
      "distinct_umis": 2,
      "reads_with_N_in_umi": 0,
      "mean_quality": 33.333333,
      "q30_fraction": 0.666667,
      "homopolymer_umi_reads": 3
    }
  }
}
```

## 7) calculate_basecalling_qc

- **Purpose:** Calculates basecalling quality metrics: no-call rate, Q30 fraction, expected errors per read and the Illumina pass-filter fraction.

- **Method:** The no-call fraction is the share of `N` bases, the Q30 fraction the share of bases with Phred Q >= 30 and the expected errors per read the mean of the sum of 10^(-Q/10) over a read's bases. Reads with a mean Q below `min_mean_quality` are counted. The pass-filter fraction comes from the CASAVA header flag (`1:Y:0:...` = filtered, `N` = passing) and is `null` when the flag is absent. Instrument-side metrics (phasing, intensity) need vendor files and are not computed. Pure Python.

- **API:** `dnakit.sequencing_qc.calculate_basecalling_qc(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input` is FASTQ text or a file path (`.gz` allowed), a list of reads ((id, sequence, quality) tuples, `{sequence, quality}` mappings, FASTQ record objects or bare sequences), or a pair `{"read1": ..., "read2": ...}` (analysed per mate). Phred qualities are required. Options: `min_mean_quality` (default 20), `phred_offset` (33 or 64, default 33), `max_reads` (default 1000000).

- **Output:** Returns a JSON-serializable envelope with `status` (`"ok"`), `module`, `operation`, `provenance` (`implementation: "native"`, `backend: "local"`), `input_summary`, `warnings` and `result`. `result` holds `reads`, `bases`, `no_call_fraction`, `reads_with_no_calls`, `mean_q`, `q30_fraction`, `expected_errors_per_read`, `min_mean_q`, `reads_below_min_mean_q`, `fraction_reads_below_min_mean_q`, `illumina_filter_flag_reads`, `pass_filter_fraction` (and `filtered_reads` when the flag is present). For a pair the result has `paired: true` and holds these fields separately under `read1` and `read2`. Invalid input or options raise `InputFormatError` or `ConfigurationError`.

- **Example code:**

```python
from dnakit.sequencing_qc import calculate_basecalling_qc

def fastq(*records):
    return "".join(f"@{name}\n{seq}\n+\n{qual}\n" for name, seq, qual in records)

reads = fastq(
    ("r1", "ACGTACGTAC", "I" * 10),
    ("r2", "ACGTACGTAC", "I" * 10),
    ("r3", "TTTTTGGGGG", "5" * 10),
    ("r4", "NNACGTACGT", "!!" + "I" * 8),
)
result = calculate_basecalling_qc(reads)
print(result)
```

- **Example output:**

```json
{
  "module": "sequencing_qc",
  "operation": "calculate_basecalling_qc",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 4,
    "text_length": 112,
    "sequence_bases": 40,
    "lengths": [10, 10, 10, 10],
    "sha256": "74e1d5397e820db4ae7e5d14647765715a6c1e7a262dff425226e4bce9a7169a"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "reads": 4,
    "bases": 40,
    "no_call_fraction": 0.05,
    "reads_with_no_calls": 1,
    "mean_q": 33.0,
    "q30_fraction": 0.7,
    "expected_errors_per_read": 0.5257,
    "min_mean_q": 20.0,
    "reads_below_min_mean_q": 0,
    "fraction_reads_below_min_mean_q": 0.0,
    "illumina_filter_flag_reads": 0,
    "pass_filter_fraction": null
  }
}
```

## 8) calculate_basecalling_qc_tenx

- **Purpose:** Calculates Cell Ranger-style basecalling metrics (Q30 in barcode, UMI and RNA read) for 10x reads.

- **Method:** Read 1 is split into barcode and UMI regions and read 2 is the RNA read; each region gets the metrics of `calculate_basecalling_qc` (no-call fraction, Q30 fraction, expected errors per read, reads below the mean-quality threshold, pass-filter fraction). The 10x chemistry is not auto-detected: `chemistry` (or the explicit lengths) decides how read 1 is split into cell barcode and UMI.

- **API:** `dnakit.sequencing_qc.calculate_basecalling_qc_tenx(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input` is a pair `{"read1": <barcode+UMI reads>, "read2": <cDNA reads>}` (FASTQ text, path or read list); a single read set is rejected and read 1 must be at least `barcode_length` + `umi_length` long. Options: `chemistry` (`v1` 14+10, `v2` 16+10, `v3` 16+12 bp barcode+UMI; default `v3`), `barcode_length`, `umi_length`, `phred_offset` (33), `max_reads` (1000000). Phred qualities are required. Extra option: `min_mean_quality` (default 20).

- **Output:** Returns a JSON-serializable envelope with `status` (`"ok"`), `module`, `operation`, `provenance` (`implementation: "native"`, `backend: "local"`), `input_summary`, `warnings` and `result`. `result` holds `chemistry`, `barcode_length`, `umi_length`, `read_pairs` and `regions` with `barcode`, `umi` and `rna_read`, each holding the fields of `calculate_basecalling_qc`. Invalid input or options raise `InputFormatError` or `ConfigurationError`.

- **Example code:**

```python
from dnakit.sequencing_qc import calculate_basecalling_qc_tenx

def fastq(*records):
    return "".join(f"@{name}\n{seq}\n+\n{qual}\n" for name, seq, qual in records)

barcode = "ACGTACGTACGTACGT"
tenx = {
    "read1": fastq(
        ("p1", barcode + "A" * 12, "I" * 28),
        ("p2", barcode + "A" * 12, "I" * 28),
        ("p3", "T" * 16 + "C" * 12, "I" * 16 + "5" * 12),
    ),
    "read2": fastq(
        ("p1", "GATTACAGATCGGAAGAGCGGTTCAGCAGG", "I" * 30),
        ("p2", "GATTACAGATCGGAAGAGCGGTTCAGCAGG", "I" * 30),
        ("p3", "C" * 30, "5" * 30),
    ),
}
result = calculate_basecalling_qc_tenx(tenx)
print(result)
```

- **Example output:**

```json
{
  "module": "sequencing_qc",
  "operation": "calculate_basecalling_qc_tenx",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 1,
    "text_length": 446,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "c2bd8879256656753e784d4d4fe4f7d2de6add20603c60c539095f494a8d364c"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "chemistry": "v3",
    "barcode_length": 16,
    "umi_length": 12,
    "read_pairs": 3,
    "regions": {
      "barcode": {
        "reads": 3,
        "bases": 48,
        "no_call_fraction": 0.0,
        "reads_with_no_calls": 0,
        "mean_q": 40.0,
        "q30_fraction": 1.0,
        "expected_errors_per_read": 0.0016,
        "min_mean_q": 20.0,
        "reads_below_min_mean_q": 0,
        "fraction_reads_below_min_mean_q": 0.0,
        "illumina_filter_flag_reads": 0,
        "pass_filter_fraction": null
      },
      "umi": {
        "reads": 3,
        "bases": 36,
        "no_call_fraction": 0.0,
        "reads_with_no_calls": 0,
        "mean_q": 33.333333,
        "q30_fraction": 0.666667,
        "expected_errors_per_read": 0.0408,
        "min_mean_q": 20.0,
        "reads_below_min_mean_q": 0,
        "fraction_reads_below_min_mean_q": 0.0,
        "illumina_filter_flag_reads": 0,
        "pass_filter_fraction": null
      },
      "rna_read": {
        "reads": 3,
        "bases": 90,
        "no_call_fraction": 0.0,
        "reads_with_no_calls": 0,
        "mean_q": 33.333333,
        "q30_fraction": 0.666667,
        "expected_errors_per_read": 0.102,
        "min_mean_q": 20.0,
        "reads_below_min_mean_q": 0,
        "fraction_reads_below_min_mean_q": 0.0,
        "illumina_filter_flag_reads": 0,
        "pass_filter_fraction": null
      }
    }
  }
}
```

## 9) calculate_contamination

- **Purpose:** Screens reads for contamination by exact k-mer containment against user-supplied reference sequences.

- **Method:** A read hits a reference when at least `min_kmer_fraction` of its canonical k-mers (strand-independent) occur in that reference. Reported are the hit counts per reference, the reads without a hit and the reads hitting several references. With `expected_reference`, the contamination fraction is the share of screened reads that hit another reference and not the expected one. This is a k-mer screen, not an alignment, and it only detects contaminants whose sequence the caller supplies. Pure Python.

- **API:** `dnakit.sequencing_qc.calculate_contamination(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input` is FASTQ text or a file path (`.gz` allowed), a list of reads ((id, sequence, quality) tuples, `{sequence, quality}` mappings, FASTQ record objects or bare sequences), or a pair `{"read1": ..., "read2": ...}` (analysed per mate). Options: `references` (required; mapping name -> sequence, or FASTA text), `kmer_size` (default 21), `min_kmer_fraction` (default 0.5), `expected_reference` (optional name), `phred_offset` (33 or 64, default 33), `max_reads` (default 1000000).

- **Output:** Returns a JSON-serializable envelope with `status` (`"ok"`), `module`, `operation`, `provenance` (`implementation: "native"`, `backend: "local"`), `input_summary`, `warnings` and `result`. `result` holds `reads_screened`, `kmer_size`, `min_kmer_fraction`, `hits` (per reference `reads` and `fraction`), `reads_without_hit`, `reads_with_multiple_hits`, `definition` and, with `expected_reference`, `expected_reference`, `contamination_fraction` and `contamination_definition`. For a pair the result has `paired: true` and holds these fields separately under `read1` and `read2`. Invalid input or options raise `InputFormatError` or `ConfigurationError`.

- **Example code:**

```python
from dnakit.sequencing_qc import calculate_contamination

host = "ACGTTGCAAGGCTTAACCGGATCGATTCGAGGCTATTGCCAGTCAAGGTACCTGATGCCGTAAGCTTGGAACCTAGGTCCA"
other = "TTGACCGGTAACCTTGGAACCGTTAGGCCATTGACCAAGGTTCCAAGGTTCCGGAATTCCGGAACCTTGGAAGGCC"
def reverse_complement(seq):
    return seq.translate(str.maketrans("ACGT", "TGCA"))[::-1]

reads = [host[:30], host[20:50], reverse_complement(host[10:40]), other[:30]]
result = calculate_contamination(
    reads,
    options={"references": {"host": host, "bug": other}, "kmer_size": 12, "expected_reference": "host"},
)
print(result)
```

- **Example output:**

```json
{
  "module": "sequencing_qc",
  "operation": "calculate_contamination",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 4,
    "text_length": 136,
    "sequence_bases": 120,
    "lengths": [30, 30, 30, 30],
    "sha256": "c2c4d13828683299b138663672c38c92457692686454ffa93ca8a805435c7942"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "reads_screened": 4,
    "kmer_size": 12,
    "min_kmer_fraction": 0.5,
    "hits": {
      "host": {
        "reads": 3,
        "fraction": 0.75
      },
      "bug": {
        "reads": 1,
        "fraction": 0.25
      }
    },
    "reads_without_hit": 0,
    "reads_with_multiple_hits": 0,
    "definition": "A read hits a reference when at least min_kmer_fraction of its canonical k-mers occur in that reference (exact k-mer containment, both strands). It is a k-mer s...",
    "expected_reference": "host",
    "contamination_fraction": 0.25,
    "contamination_definition": "Fraction of screened reads that hit at least one non-expected reference and not the expected one."
  }
}
```

## 10) calculate_contamination_tenx

- **Purpose:** Screens the cDNA reads (read 2) of 10x data for contamination by k-mer containment against supplied references.

- **Method:** Read 2 is screened exactly as in `calculate_contamination` (canonical k-mers, `min_kmer_fraction` rule, optional `expected_reference`); read 1 only supplies the layout and the read count. The 10x chemistry is not auto-detected: `chemistry` (or the explicit lengths) decides how read 1 is split into cell barcode and UMI.

- **API:** `dnakit.sequencing_qc.calculate_contamination_tenx(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input` is a pair `{"read1": <barcode+UMI reads>, "read2": <cDNA reads>}` (FASTQ text, path or read list); a single read set is rejected and read 1 must be at least `barcode_length` + `umi_length` long. Options: `chemistry` (`v1` 14+10, `v2` 16+10, `v3` 16+12 bp barcode+UMI; default `v3`), `barcode_length`, `umi_length`, `phred_offset` (33), `max_reads` (1000000). Extra options as `calculate_contamination`: `references` (required), `kmer_size` (21), `min_kmer_fraction` (0.5), `expected_reference`.

- **Output:** Returns a JSON-serializable envelope with `status` (`"ok"`), `module`, `operation`, `provenance` (`implementation: "native"`, `backend: "local"`), `input_summary`, `warnings` and `result`. `result` holds `chemistry`, `barcode_length`, `umi_length`, `read_pairs`, `region` (`rna_read (read 2)`) and the fields of `calculate_contamination` (`reads_screened`, `hits`, `reads_without_hit`, `reads_with_multiple_hits`, ...). Invalid input or options raise `InputFormatError` or `ConfigurationError`.

- **Example code:**

```python
from dnakit.sequencing_qc import calculate_contamination_tenx

def fastq(*records):
    return "".join(f"@{name}\n{seq}\n+\n{qual}\n" for name, seq, qual in records)

barcode = "ACGTACGTACGTACGT"
tenx = {
    "read1": fastq(
        ("p1", barcode + "A" * 12, "I" * 28),
        ("p2", barcode + "A" * 12, "I" * 28),
        ("p3", "T" * 16 + "C" * 12, "I" * 16 + "5" * 12),
    ),
    "read2": fastq(
        ("p1", "GATTACAGATCGGAAGAGCGGTTCAGCAGG", "I" * 30),
        ("p2", "GATTACAGATCGGAAGAGCGGTTCAGCAGG", "I" * 30),
        ("p3", "C" * 30, "5" * 30),
    ),
}
host = "GATTACAGATCGGAAGAGCGGTTCAGCAGGTTACCAGTCAGGCCTAAGGTCCATTGGCC"
result = calculate_contamination_tenx(tenx, options={"references": {"host": host}, "kmer_size": 12})
print(result)
```

- **Example output:**

```json
{
  "module": "sequencing_qc",
  "operation": "calculate_contamination_tenx",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 1,
    "text_length": 446,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "c2bd8879256656753e784d4d4fe4f7d2de6add20603c60c539095f494a8d364c"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "chemistry": "v3",
    "barcode_length": 16,
    "umi_length": 12,
    "read_pairs": 3,
    "region": "rna_read (read 2)",
    "reads_screened": 3,
    "kmer_size": 12,
    "min_kmer_fraction": 0.5,
    "hits": {
      "host": {
        "reads": 2,
        "fraction": 0.666667
      }
    },
    "reads_without_hit": 1,
    "reads_with_multiple_hits": 0,
    "definition": "A read hits a reference when at least min_kmer_fraction of its canonical k-mers occur in that reference (exact k-mer containment, both strands). It is a k-mer s..."
  }
}
```

## 11) calculate_duplicate_rate

- **Purpose:** Calculates the exact-sequence duplicate rate of reads: (reads - distinct sequences) / reads.

- **Method:** Reads are compared as exact sequences (optionally only the first `prefix_length` bases, optionally treating a read and its reverse complement as the same). For a pair the key is the concatenation of both mates. Read positions are unknown, so this is sequence duplication, not mapping duplication as marked by Picard. Pure Python.

- **API:** `dnakit.sequencing_qc.calculate_duplicate_rate(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input` is FASTQ text or a file path (`.gz` allowed), a list of reads ((id, sequence, quality) tuples, `{sequence, quality}` mappings, FASTQ record objects or bare sequences), or a pair `{"read1": ..., "read2": ...}` (one key per pair). Options: `prefix_length` (compare only the first N bases; default 0 = whole read), `reverse_complement_aware` (default `False`), `phred_offset` (33 or 64, default 33), `max_reads` (default 1000000).

- **Output:** Returns a JSON-serializable envelope with `status` (`"ok"`), `module`, `operation`, `provenance` (`implementation: "native"`, `backend: "local"`), `input_summary`, `warnings` and `result`. `result` holds `paired`, `reads`, `distinct_sequences`, `duplicate_reads`, `duplicate_rate`, `unique_fraction`, `most_common` (sequence and count) and a `definition` string. Invalid input or options raise `InputFormatError` or `ConfigurationError`.

- **Example code:**

```python
from dnakit.sequencing_qc import calculate_duplicate_rate

def fastq(*records):
    return "".join(f"@{name}\n{seq}\n+\n{qual}\n" for name, seq, qual in records)

reads = fastq(
    ("r1", "ACGTACGTAC", "I" * 10),
    ("r2", "ACGTACGTAC", "I" * 10),
    ("r3", "TTTTTGGGGG", "5" * 10),
    ("r4", "NNACGTACGT", "!!" + "I" * 8),
)
result = calculate_duplicate_rate(reads)
print(result)
```

- **Example output:**

```json
{
  "module": "sequencing_qc",
  "operation": "calculate_duplicate_rate",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 4,
    "text_length": 112,
    "sequence_bases": 40,
    "lengths": [10, 10, 10, 10],
    "sha256": "74e1d5397e820db4ae7e5d14647765715a6c1e7a262dff425226e4bce9a7169a"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "paired": false,
    "reads": 4,
    "distinct_sequences": 3,
    "duplicate_reads": 1,
    "duplicate_rate": 0.25,
    "unique_fraction": 0.75,
    "most_common": [
      {
        "sequence": "ACGTACGTAC",
        "count": 2
      },
      {
        "sequence": "NNACGTACGT",
        "count": 1
      },
      "..."
    ],
    "definition": "duplicate_rate = (reads - distinct sequences) / reads; sequences are compared exactly (read pairs by the concatenation of both mates)."
  }
}
```

## 12) calculate_duplication

- **Purpose:** Calculates the FastQC-style sequence duplication level histogram of reads.

- **Method:** The level of a sequence is the number of copies of that exact sequence, binned as 1-9, >10 (10-49), >50, >100, >500, >1k, >5k and >10k+ as in FastQC. Each level is reported as a percentage of the total reads and of the distinct sequences, together with the percentage of reads that remain after deduplication. No FastQC-style extrapolation to unseen reads is applied. Pure Python.

- **API:** `dnakit.sequencing_qc.calculate_duplication(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input` is FASTQ text or a file path (`.gz` allowed), a list of reads ((id, sequence, quality) tuples, `{sequence, quality}` mappings, FASTQ record objects or bare sequences), or a pair `{"read1": ..., "read2": ...}` (one key per pair). Options: `prefix_length` (default 0), `reverse_complement_aware` (default `False`), `phred_offset` (33 or 64, default 33), `max_reads` (default 1000000).

- **Output:** Returns a JSON-serializable envelope with `status` (`"ok"`), `module`, `operation`, `provenance` (`implementation: "native"`, `backend: "local"`), `input_summary`, `warnings` and `result`. `result` holds `paired`, `reads`, `distinct_sequences`, `percent_remaining_after_deduplication`, `levels` (16 rows with `level`, `distinct_sequences`, `percent_of_total_reads`, `percent_of_deduplicated`) and a `definition` string. Invalid input or options raise `InputFormatError` or `ConfigurationError`.

- **Example code:**

```python
from dnakit.sequencing_qc import calculate_duplication

reads = ["AAAA"] * 3 + ["CCCC", "GGGG"]
result = calculate_duplication(reads)
print(result)
```

- **Example output:**

```json
{
  "module": "sequencing_qc",
  "operation": "calculate_duplication",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 5,
    "text_length": 40,
    "sequence_bases": 20,
    "lengths": [4, 4, 4, 4, 4],
    "sha256": "9a2abe2e853265a2f1ad0d543707ac64ac86a065f17d998ffdb72db0f8f0dc35"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "paired": false,
    "reads": 5,
    "distinct_sequences": 3,
    "percent_remaining_after_deduplication": 60.0,
    "levels": [
      {
        "level": "1",
        "distinct_sequences": 2,
        "percent_of_total_reads": 40.0,
        "percent_of_deduplicated": 66.666667
      },
      {
        "level": "2",
        "distinct_sequences": 0,
        "percent_of_total_reads": 0.0,
        "percent_of_deduplicated": 0.0
      },
      "..."
    ],
    "definition": "Level = number of copies of an identical sequence; bins 1-9, >10 (10-49), >50, >100, >500, >1k, >5k, >10k+ as in FastQC. No extrapolation beyond the supplied re..."
  }
}
```

## 13) calculate_duplication_tenx

- **Purpose:** Calculates the duplication of 10x reads, keyed on cell barcode + UMI + the first bases of read 2.

- **Method:** Two reads are duplicates when their cell barcode, UMI and the first `cdna_prefix_length` bases of read 2 are identical. This is a heuristic proxy for sequencing saturation: Cell Ranger uses (barcode, UMI, gene) from an alignment, which is unavailable here. The levels are binned as in `calculate_duplication`. The 10x chemistry is not auto-detected: `chemistry` (or the explicit lengths) decides how read 1 is split into cell barcode and UMI.

- **API:** `dnakit.sequencing_qc.calculate_duplication_tenx(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input` is a pair `{"read1": <barcode+UMI reads>, "read2": <cDNA reads>}` (FASTQ text, path or read list); a single read set is rejected and read 1 must be at least `barcode_length` + `umi_length` long. Options: `chemistry` (`v1` 14+10, `v2` 16+10, `v3` 16+12 bp barcode+UMI; default `v3`), `barcode_length`, `umi_length`, `phred_offset` (33), `max_reads` (1000000). Extra option: `cdna_prefix_length` (default 25).

- **Output:** Returns a JSON-serializable envelope with `status` (`"ok"`), `module`, `operation`, `provenance` (`implementation: "native"`, `backend: "local"`), `input_summary`, `warnings` and `result`. `result` holds `chemistry`, `barcode_length`, `umi_length`, `read_pairs`, `reads`, `distinct_sequences`, `percent_remaining_after_deduplication`, `levels`, `definition`, `duplicate_rate` and `sequencing_saturation_proxy`. Invalid input or options raise `InputFormatError` or `ConfigurationError`.

- **Example code:**

```python
from dnakit.sequencing_qc import calculate_duplication_tenx

def fastq(*records):
    return "".join(f"@{name}\n{seq}\n+\n{qual}\n" for name, seq, qual in records)

barcode = "ACGTACGTACGTACGT"
tenx = {
    "read1": fastq(
        ("p1", barcode + "A" * 12, "I" * 28),
        ("p2", barcode + "A" * 12, "I" * 28),
        ("p3", "T" * 16 + "C" * 12, "I" * 16 + "5" * 12),
    ),
    "read2": fastq(
        ("p1", "GATTACAGATCGGAAGAGCGGTTCAGCAGG", "I" * 30),
        ("p2", "GATTACAGATCGGAAGAGCGGTTCAGCAGG", "I" * 30),
        ("p3", "C" * 30, "5" * 30),
    ),
}
result = calculate_duplication_tenx(tenx)
print(result)
```

- **Example output:**

```json
{
  "module": "sequencing_qc",
  "operation": "calculate_duplication_tenx",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 1,
    "text_length": 446,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "c2bd8879256656753e784d4d4fe4f7d2de6add20603c60c539095f494a8d364c"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "reads": 3,
    "distinct_sequences": 2,
    "percent_remaining_after_deduplication": 66.666667,
    "levels": [
      {
        "level": "1",
        "distinct_sequences": 1,
        "percent_of_total_reads": 33.333333,
        "percent_of_deduplicated": 50.0
      },
      {
        "level": "2",
        "distinct_sequences": 1,
        "percent_of_total_reads": 66.666667,
        "percent_of_deduplicated": 50.0
      },
      "..."
    ],
    "definition": "Reads are duplicates when cell barcode, UMI and the first 25 bases of read 2 are identical. This is a proxy: Cell Ranger's sequencing saturation uses (barcode, ...",
    "chemistry": "v3",
    "barcode_length": 16,
    "umi_length": 12,
    "read_pairs": 3,
    "duplicate_rate": 0.333333,
    "sequencing_saturation_proxy": 0.333333
  }
}
```

## 14) calculate_error_profile

- **Purpose:** Calculates the per-position expected error profile of reads from their Phred scores.

- **Method:** The expected error rate at a position is the mean of 10^(-Q/10) over the reads that cover it; the result also gives the cumulative expected errors along the read, the expected errors per read and the `N` fraction per position. This is the error the quality scores claim, not an observed mismatch rate (use `estimate_read_error_model` with a reference for that). Pure Python.

- **API:** `dnakit.sequencing_qc.calculate_error_profile(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input` is FASTQ text or a file path (`.gz` allowed), a list of reads ((id, sequence, quality) tuples, `{sequence, quality}` mappings, FASTQ record objects or bare sequences), or a pair `{"read1": ..., "read2": ...}` (analysed per mate). Phred qualities are required. Options: `phred_offset` (33 or 64, default 33), `max_reads` (default 1000000).

- **Output:** Returns a JSON-serializable envelope with `status` (`"ok"`), `module`, `operation`, `provenance` (`implementation: "native"`, `backend: "local"`), `input_summary`, `warnings` and `result`. `result` holds `source` (`phred_quality`), `overall_expected_error_rate`, `overall_q_equivalent`, `expected_errors_per_read`, `per_position` (`position`, `bases`, `expected_error_rate`, `cumulative_expected_errors`, `n_fraction`) and a `note`. For a pair the result has `paired: true` and holds these fields separately under `read1` and `read2`. Invalid input or options raise `InputFormatError` or `ConfigurationError`.

- **Example code:**

```python
from dnakit.sequencing_qc import calculate_error_profile

def fastq(*records):
    return "".join(f"@{name}\n{seq}\n+\n{qual}\n" for name, seq, qual in records)

reads = fastq(
    ("r1", "ACGTACGTAC", "I" * 10),
    ("r2", "ACGTACGTAC", "I" * 10),
    ("r3", "TTTTTGGGGG", "5" * 10),
    ("r4", "NNACGTACGT", "!!" + "I" * 8),
)
result = calculate_error_profile(reads)
print(result)
```

- **Example output:**

```json
{
  "module": "sequencing_qc",
  "operation": "calculate_error_profile",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 4,
    "text_length": 112,
    "sequence_bases": 40,
    "lengths": [10, 10, 10, 10],
    "sha256": "74e1d5397e820db4ae7e5d14647765715a6c1e7a262dff425226e4bce9a7169a"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "source": "phred_quality",
    "overall_expected_error_rate": 0.05256999999999999,
    "overall_q_equivalent": 12.79262,
    "expected_errors_per_read": 0.5257,
    "per_position": [
      {
        "position": 1,
        "bases": 4,
        "expected_error_rate": 0.25255,
        "cumulative_expected_errors": 0.25255,
        "n_fraction": 0.25
      },
      {
        "position": 2,
        "bases": 4,
        "expected_error_rate": 0.25255,
        "cumulative_expected_errors": 0.5051,
        "n_fraction": 0.25
      },
      "..."
    ],
    "note": "Expected error = 10^(-Q/10) from the reported Phred scores; not an observed mismatch rate."
  }
}
```

## 15) calculate_error_profile_tenx

- **Purpose:** Calculates Phred-derived expected error profiles for the barcode, UMI and RNA regions of 10x reads.

- **Method:** Read 1 is split into barcode and UMI regions and read 2 is the RNA read; each region gets the result of `calculate_error_profile`, i.e. the expected error rate 10^(-Q/10) averaged per position. The 10x chemistry is not auto-detected: `chemistry` (or the explicit lengths) decides how read 1 is split into cell barcode and UMI.

- **API:** `dnakit.sequencing_qc.calculate_error_profile_tenx(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input` is a pair `{"read1": <barcode+UMI reads>, "read2": <cDNA reads>}` (FASTQ text, path or read list); a single read set is rejected and read 1 must be at least `barcode_length` + `umi_length` long. Options: `chemistry` (`v1` 14+10, `v2` 16+10, `v3` 16+12 bp barcode+UMI; default `v3`), `barcode_length`, `umi_length`, `phred_offset` (33), `max_reads` (1000000). Phred qualities are required.

- **Output:** Returns a JSON-serializable envelope with `status` (`"ok"`), `module`, `operation`, `provenance` (`implementation: "native"`, `backend: "local"`), `input_summary`, `warnings` and `result`. `result` holds `chemistry`, `barcode_length`, `umi_length`, `read_pairs` and `regions` with `barcode`, `umi` and `rna_read`, each holding the fields of `calculate_error_profile`. Invalid input or options raise `InputFormatError` or `ConfigurationError`.

- **Example code:**

```python
from dnakit.sequencing_qc import calculate_error_profile_tenx

def fastq(*records):
    return "".join(f"@{name}\n{seq}\n+\n{qual}\n" for name, seq, qual in records)

barcode = "ACGTACGTACGTACGT"
tenx = {
    "read1": fastq(
        ("p1", barcode + "A" * 12, "I" * 28),
        ("p2", barcode + "A" * 12, "I" * 28),
        ("p3", "T" * 16 + "C" * 12, "I" * 16 + "5" * 12),
    ),
    "read2": fastq(
        ("p1", "GATTACAGATCGGAAGAGCGGTTCAGCAGG", "I" * 30),
        ("p2", "GATTACAGATCGGAAGAGCGGTTCAGCAGG", "I" * 30),
        ("p3", "C" * 30, "5" * 30),
    ),
}
result = calculate_error_profile_tenx(tenx)
print(result)
```

- **Example output:**

```json
{
  "module": "sequencing_qc",
  "operation": "calculate_error_profile_tenx",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 1,
    "text_length": 446,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "c2bd8879256656753e784d4d4fe4f7d2de6add20603c60c539095f494a8d364c"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "chemistry": "v3",
    "barcode_length": 16,
    "umi_length": 12,
    "read_pairs": 3,
    "regions": {
      "barcode": {
        "source": "phred_quality",
        "overall_expected_error_rate": 0.0001,
        "overall_q_equivalent": 40.0,
        "expected_errors_per_read": 0.0016,
        "per_position": [
          {
            "...": "5 keys"
          },
          {
            "...": "5 keys"
          },
          "..."
        ],
        "note": "Expected error = 10^(-Q/10) from the reported Phred scores; not an observed mismatch rate."
      },
      "umi": {
        "source": "phred_quality",
        "overall_expected_error_rate": 0.0034000000000000002,
        "overall_q_equivalent": 24.685211,
        "expected_errors_per_read": 0.0408,
        "per_position": [
          {
            "...": "5 keys"
          },
          {
            "...": "5 keys"
          },
          "..."
        ],
        "note": "Expected error = 10^(-Q/10) from the reported Phred scores; not an observed mismatch rate."
      },
      "rna_read": {
        "source": "phred_quality",
        "overall_expected_error_rate": 0.0034000000000000007,
        "overall_q_equivalent": 24.685211,
        "expected_errors_per_read": 0.10200000000000002,
        "per_position": [
          {
            "...": "5 keys"
          },
          {
            "...": "5 keys"
          },
          "..."
        ],
        "note": "Expected error = 10^(-Q/10) from the reported Phred scores; not an observed mismatch rate."
      }
    }
  }
}
```

## 16) calculate_gc_bias

- **Purpose:** Calculates the GC content distribution of reads and its deviation from a fitted normal distribution (FastQC convention).

- **Method:** The GC percentage of every read (G+C over A+C+G+T) is histogrammed in 1% bins and compared with a normal distribution having the observed mean and standard deviation. `deviation_from_normal_percent` is 100 x the sum over bins of |observed - theoretical count| divided by the number of reads (`null` when the standard deviation is 0). GC% by read position is also given. With `reference_gc` the shift of the mean from the expected GC is reported. This is a read-composition profile, not Picard's coverage-versus-GC bias, which needs alignments. Pure Python.

- **API:** `dnakit.sequencing_qc.calculate_gc_bias(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input` is FASTQ text or a file path (`.gz` allowed), a list of reads ((id, sequence, quality) tuples, `{sequence, quality}` mappings, FASTQ record objects or bare sequences), or a pair `{"read1": ..., "read2": ...}` (analysed per mate). Options: `reference_gc` (optional expected GC fraction in [0, 1]), `phred_offset` (33 or 64, default 33), `max_reads` (default 1000000).

- **Output:** Returns a JSON-serializable envelope with `status` (`"ok"`), `module`, `operation`, `provenance` (`implementation: "native"`, `backend: "local"`), `input_summary`, `warnings` and `result`. `result` holds `reads_with_gc`, `mean_gc_percent`, `stdev_gc_percent`, `mode_gc_percent`, `histogram` (GC% -> read count), `theoretical_normal_counts`, `deviation_from_normal_percent`, `gc_percent_by_position`, `definition` and, with `reference_gc`, `reference_gc_percent` and `gc_shift_vs_reference_percent`. For a pair the result has `paired: true` and holds these fields separately under `read1` and `read2`. Invalid input or options raise `InputFormatError` or `ConfigurationError`.

- **Example code:**

```python
from dnakit.sequencing_qc import calculate_gc_bias

reads = ["GGGG", "AAAA", "GCAT", "GGAT"]
result = calculate_gc_bias(reads, options={"reference_gc": 0.4})
print(result)
```

- **Example output:**

```json
{
  "module": "sequencing_qc",
  "operation": "calculate_gc_bias",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 4,
    "text_length": 32,
    "sequence_bases": 16,
    "lengths": [4, 4, 4, 4],
    "sha256": "2a486830d70f56db6c17aebf31a0bb6fc9ec7c1df90a1cb1dfbeced28e260d66"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "reads_with_gc": 4,
    "mean_gc_percent": 50.0,
    "stdev_gc_percent": 40.824829,
    "mode_gc_percent": 50,
    "histogram": {
      "0": 1,
      "50": 2,
      "100": 1
    },
    "theoretical_normal_counts": [
      0.0185,
      0.019,
      "..."
    ],
    "deviation_from_normal_percent": 174.590226,
    "gc_percent_by_position": [
      75.0,
      75.0,
      "..."
    ],
    "definition": "Per-read GC% (G+C over A+C+G+T) histogrammed in 1% bins; the theoretical curve is a normal distribution with the observed mean and standard deviation (FastQC co...",
    "reference_gc_percent": 40.0,
    "gc_shift_vs_reference_percent": 10.0
  }
}
```

## 17) calculate_gc_bias_tenx

- **Purpose:** Calculates the GC content distribution of 10x cDNA reads (read 2).

- **Method:** Read 2 is profiled exactly as in `calculate_gc_bias` (1% bins, fitted normal curve, deviation, GC% by position). The 10x chemistry is not auto-detected: `chemistry` (or the explicit lengths) decides how read 1 is split into cell barcode and UMI.

- **API:** `dnakit.sequencing_qc.calculate_gc_bias_tenx(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input` is a pair `{"read1": <barcode+UMI reads>, "read2": <cDNA reads>}` (FASTQ text, path or read list); a single read set is rejected and read 1 must be at least `barcode_length` + `umi_length` long. Options: `chemistry` (`v1` 14+10, `v2` 16+10, `v3` 16+12 bp barcode+UMI; default `v3`), `barcode_length`, `umi_length`, `phred_offset` (33), `max_reads` (1000000). Extra option: `reference_gc` (optional GC fraction in [0, 1]).

- **Output:** Returns a JSON-serializable envelope with `status` (`"ok"`), `module`, `operation`, `provenance` (`implementation: "native"`, `backend: "local"`), `input_summary`, `warnings` and `result`. `result` holds `chemistry`, `barcode_length`, `umi_length`, `read_pairs`, `region` (`rna_read (read 2)`) and the fields of `calculate_gc_bias`. Invalid input or options raise `InputFormatError` or `ConfigurationError`.

- **Example code:**

```python
from dnakit.sequencing_qc import calculate_gc_bias_tenx

def fastq(*records):
    return "".join(f"@{name}\n{seq}\n+\n{qual}\n" for name, seq, qual in records)

barcode = "ACGTACGTACGTACGT"
tenx = {
    "read1": fastq(
        ("p1", barcode + "A" * 12, "I" * 28),
        ("p2", barcode + "A" * 12, "I" * 28),
        ("p3", "T" * 16 + "C" * 12, "I" * 16 + "5" * 12),
    ),
    "read2": fastq(
        ("p1", "GATTACAGATCGGAAGAGCGGTTCAGCAGG", "I" * 30),
        ("p2", "GATTACAGATCGGAAGAGCGGTTCAGCAGG", "I" * 30),
        ("p3", "C" * 30, "5" * 30),
    ),
}
result = calculate_gc_bias_tenx(tenx, options={"reference_gc": 0.4})
print(result)
```

- **Example output:**

```json
{
  "module": "sequencing_qc",
  "operation": "calculate_gc_bias_tenx",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 1,
    "text_length": 446,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "c2bd8879256656753e784d4d4fe4f7d2de6add20603c60c539095f494a8d364c"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "chemistry": "v3",
    "barcode_length": 16,
    "umi_length": 12,
    "read_pairs": 3,
    "region": "rna_read (read 2)",
    "reads_with_gc": 3,
    "mean_gc_percent": 68.888889,
    "stdev_gc_percent": 26.943013,
    "mode_gc_percent": 53,
    "histogram": {
      "53": 2,
      "100": 1
    },
    "theoretical_normal_counts": [
      0.0017,
      0.0019,
      "..."
    ],
    "deviation_from_normal_percent": 183.455679,
    "gc_percent_by_position": [
      100.0,
      33.333333,
      "..."
    ],
    "definition": "Per-read GC% (G+C over A+C+G+T) histogrammed in 1% bins; the theoretical curve is a normal distribution with the observed mean and standard deviation (FastQC co...",
    "reference_gc_percent": 40.0,
    "gc_shift_vs_reference_percent": 28.888889
  }
}
```

## 18) calculate_index_hopping

- **Purpose:** Calculates the index-hopping rate of dual-indexed reads from the sample indexes in their headers.

- **Method:** The `i7+i5` index is read from the end of the Illumina header comment (for example `1:N:0:ACGTAC+TTGGCA`). A read is hopped when its i7 and i5 each occur in the sample sheet but not as a pair; rate = hopped / (hopped + reads with an expected pair). Indexes must match exactly (no mismatch tolerance). Reads whose indexes are not in the sheet are reported separately as unknown and are not part of the rate. Pure Python.

- **API:** `dnakit.sequencing_qc.calculate_index_hopping(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input` is FASTQ text or a file path (`.gz` allowed), a list of reads ((id, sequence, quality) tuples, `{sequence, quality}` mappings, FASTQ record objects or bare sequences) whose header comment ends in `i7+i5`. A pair mapping is not accepted (`INVALID_INPUT`); for paired data pass one mate (both mates carry the same indexes). Options: `samples` (required; mapping sample -> `i7+i5`), `phred_offset` (33 or 64, default 33), `max_reads` (default 1000000).

- **Output:** Returns a JSON-serializable envelope with `status` (`"ok"`), `module`, `operation`, `provenance` (`implementation: "native"`, `backend: "local"`), `input_summary`, `warnings` and `result`. `result` holds `reads`, `expected_combination_reads`, `index_hopped_reads`, `unknown_index_reads`, `index_hopping_rate`, `unknown_index_fraction`, `hopped_combinations` (`index`, `reads`), `reads_per_expected_sample` and a `definition` string. Invalid input or options raise `InputFormatError` or `ConfigurationError`.

- **Example code:**

```python
from dnakit.sequencing_qc import calculate_index_hopping

records = (
    [("s1", "A", "I", "1:N:0:AAAA+CCCC")] * 6
    + [("s2", "A", "I", "1:N:0:GGGG+TTTT")] * 2
    + [("hop", "A", "I", "1:N:0:AAAA+TTTT")]
    + [("unk", "A", "I", "1:N:0:CCCC+AAAA")]
)
reads = "".join(f"@{name} {comment}\n{seq}\n+\n{qual}\n" for name, seq, qual, comment in records)
result = calculate_index_hopping(reads, options={"samples": {"S1": "AAAA+CCCC", "S2": "GGGG+TTTT"}})
print(result)
```

- **Example output:**

```json
{
  "module": "sequencing_qc",
  "operation": "calculate_index_hopping",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 10,
    "text_length": 262,
    "sequence_bases": 10,
    "lengths": [1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
    "sha256": "3be23d01c8b0b3436149a94c2fe9a3fb2a272880ae45bc5fdf7650fe52f6b2ee"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "reads": 10,
    "expected_combination_reads": 8,
    "index_hopped_reads": 1,
    "unknown_index_reads": 1,
    "index_hopping_rate": 0.111111,
    "unknown_index_fraction": 0.1,
    "hopped_combinations": [
      {
        "index": "AAAA+TTTT",
        "reads": 1
      }
    ],
    "reads_per_expected_sample": {
      "S1": 6,
      "S2": 2
    },
    "definition": "A read is index-hopped when its i7 and i5 are each used in the sample sheet but the pair is not. rate = hopped / (hopped + expected-pair reads). Reads whose ind..."
  }
}
```

## 19) calculate_index_hopping_tenx

- **Purpose:** Calculates the index-hopping rate of dual-indexed 10x libraries from the sample indexes in read headers.

- **Method:** Same computation as `calculate_index_hopping` (exact `i7+i5` matching, rate = hopped / (hopped + expected-pair reads)). It applies to 10x dual-index kits; single-index kits carry one index per sample and no hopping signal, and are rejected with an error.

- **API:** `dnakit.sequencing_qc.calculate_index_hopping_tenx(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input` is FASTQ text or a file path (`.gz` allowed), a list of reads ((id, sequence, quality) tuples, `{sequence, quality}` mappings, FASTQ record objects or bare sequences) whose header comment ends in `i7+i5`. A pair mapping is not accepted (`INVALID_INPUT`); for paired data pass one mate. Options as `calculate_index_hopping`: `samples` (required), `phred_offset` (33 or 64, default 33), `max_reads` (default 1000000).

- **Output:** Returns a JSON-serializable envelope with `status` (`"ok"`), `module`, `operation`, `provenance` (`implementation: "native"`, `backend: "local"`), `input_summary`, `warnings` and `result`. `result` holds the same fields as `calculate_index_hopping` (`reads`, `expected_combination_reads`, `index_hopped_reads`, `unknown_index_reads`, `index_hopping_rate`, `unknown_index_fraction`, `hopped_combinations`, `reads_per_expected_sample`, `definition`). Invalid input or options raise `InputFormatError` or `ConfigurationError`.

- **Example code:**

```python
from dnakit.sequencing_qc import calculate_index_hopping_tenx

records = (
    [("s1", "A", "I", "1:N:0:AAAA+CCCC")] * 6
    + [("s2", "A", "I", "1:N:0:GGGG+TTTT")] * 2
    + [("hop", "A", "I", "1:N:0:AAAA+TTTT")]
    + [("unk", "A", "I", "1:N:0:CCCC+AAAA")]
)
reads = "".join(f"@{name} {comment}\n{seq}\n+\n{qual}\n" for name, seq, qual, comment in records)
result = calculate_index_hopping_tenx(reads, options={"samples": {"S1": "AAAA+CCCC", "S2": "GGGG+TTTT"}})
print(result)
```

- **Example output:**

```json
{
  "module": "sequencing_qc",
  "operation": "calculate_index_hopping_tenx",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 10,
    "text_length": 262,
    "sequence_bases": 10,
    "lengths": [1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
    "sha256": "3be23d01c8b0b3436149a94c2fe9a3fb2a272880ae45bc5fdf7650fe52f6b2ee"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "reads": 10,
    "expected_combination_reads": 8,
    "index_hopped_reads": 1,
    "unknown_index_reads": 1,
    "index_hopping_rate": 0.111111,
    "unknown_index_fraction": 0.1,
    "hopped_combinations": [
      {
        "index": "AAAA+TTTT",
        "reads": 1
      }
    ],
    "reads_per_expected_sample": {
      "S1": 6,
      "S2": 2
    },
    "definition": "A read is index-hopped when its i7 and i5 are each used in the sample sheet but the pair is not. rate = hopped / (hopped + expected-pair reads). Reads whose ind..."
  }
}
```

## 20) calculate_insert_size

- **Purpose:** Calculates insert (fragment) size statistics from read-pair overlap, SAM text, a 10x fragments table or a list of sizes.

- **Method:** Depending on the input the sizes come from the overlap of read 1 with the reverse complement of read 2 (pairs without a detectable overlap, i.e. an insert longer than both reads combined, are excluded and counted), from the absolute template length (TLEN) of mapped first-in-pair SAM reads (proper pairs only by default), from the fragments table (end - start, weighted by count) or directly from the numbers. Sizes above `max_insert_size` are excluded and counted. Statistics are the count, range, mean, median, standard deviation, median absolute deviation, mode and percentiles. Pure Python.

- **API:** `dnakit.sequencing_qc.calculate_insert_size(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input` is a read pair `{"read1": ..., "read2": ...}` (size = overlap-based fragment length, i.e. length of read 1 + read 2 minus their overlap, where read 1 overlaps the reverse complement of read 2), SAM text or path (|TLEN| of mapped first-in-pair reads), a fragments/BED-like TSV (`chrom`, `start`, `end` and optionally `barcode`, `count`; size = end - start, weighted by count) or a list of numbers. Options: `min_overlap` (default 10), `max_mismatch_rate` (default 0.1), `max_pairs` (pairs sized by overlap, default 200000), `proper_pairs_only` (SAM, default `True`), `max_insert_size` (default 100000), `phred_offset`, `max_reads`.

- **Output:** Returns a JSON-serializable envelope with `status` (`"ok"`), `module`, `operation`, `provenance` (`implementation: "native"`, `backend: "local"`), `input_summary`, `warnings` and `result`. `result` holds `source` (`read-pair overlap`, `SAM TLEN`, `fragments (end - start)` or `numbers`), `count`, `excluded_above_max`, `min`, `max`, `mean`, `median`, `stdev`, `median_absolute_deviation`, `mode`, `percentiles` (`p5`, `p25`, `p75`, `p95`) and source-specific counts (`pairs`, `pairs_without_overlap`, `non_proper_pairs_skipped`, ...). Invalid input or options raise `InputFormatError` or `ConfigurationError`.

- **Example code:**

```python
from dnakit.sequencing_qc import calculate_insert_size

fragment = "GATTACAGCTAGGCTTACGGATCCATGCAA"

def reverse_complement(seq):
    return seq.translate(str.maketrans("ACGT", "TGCA"))[::-1]

pairs = {"read1": [fragment[:20]], "read2": [reverse_complement(fragment[10:])]}
result = calculate_insert_size(pairs, options={"min_overlap": 10, "max_mismatch_rate": 0.0})
print(result)
```

- **Example output:**

```json
{
  "module": "sequencing_qc",
  "operation": "calculate_insert_size",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 1,
    "text_length": 70,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "52677fd537fe6e846e979fc7d0a79c8e4dfe046926d05023ec12df1cfd1de470"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "source": "read-pair overlap",
    "pairs": 1,
    "pairs_without_overlap": 0,
    "min_overlap": 10,
    "max_mismatch_rate": 0.0,
    "count": 1,
    "excluded_above_max": 0,
    "min": 30,
    "max": 30,
    "mean": 30.0,
    "median": 30.0,
    "stdev": null,
    "median_absolute_deviation": 0.0,
    "mode": 30,
    "percentiles": {
      "p5": 30.0,
      "p25": 30.0,
      "p75": 30.0,
      "p95": 30.0
    }
  }
}
```

## 21) calculate_insert_size_distribution

- **Purpose:** Calculates the insert (fragment) size histogram with cumulative fractions.

- **Method:** Sizes are obtained as in `calculate_insert_size` and then binned into windows of `bin_width` (default 1, widened automatically for wide ranges). Each bin gives its count, fraction of all sizes and cumulative fraction; `end` is the last size inside the bin.

- **API:** `dnakit.sequencing_qc.calculate_insert_size_distribution(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input` and options as `calculate_insert_size`, plus `bin_width` (default 1). `input` is a read pair `{"read1": ..., "read2": ...}` (size = overlap-based fragment length, i.e. length of read 1 + read 2 minus their overlap, where read 1 overlaps the reverse complement of read 2), SAM text or path (|TLEN| of mapped first-in-pair reads), a fragments/BED-like TSV (`chrom`, `start`, `end` and optionally `barcode`, `count`; size = end - start, weighted by count) or a list of numbers. Options: `min_overlap` (default 10), `max_mismatch_rate` (default 0.1), `max_pairs` (pairs sized by overlap, default 200000), `proper_pairs_only` (SAM, default `True`), `max_insert_size` (default 100000), `phred_offset`, `max_reads`.

- **Output:** Returns a JSON-serializable envelope with `status` (`"ok"`), `module`, `operation`, `provenance` (`implementation: "native"`, `backend: "local"`), `input_summary`, `warnings` and `result`. `result` holds the summary statistics of `calculate_insert_size` (`count`, `excluded_above_max`, `min`, `max`, `mean`, `median`, `stdev`, `median_absolute_deviation`, `mode`, `percentiles` (`p5`, `p25`, `p75`, `p95`) and source-specific counts (`pairs`, `pairs_without_overlap`, `non_proper_pairs_skipped`, ...)), `bin_width` and `histogram` (per bin `start`, `end`, `count`, `fraction`, `cumulative_fraction`). Invalid input or options raise `InputFormatError` or `ConfigurationError`.

- **Example code:**

```python
from dnakit.sequencing_qc import calculate_insert_size_distribution

result = calculate_insert_size_distribution([100, 100, 200, 300], options={"bin_width": 100})
print(result)
```

- **Example output:**

```json
{
  "module": "sequencing_qc",
  "operation": "calculate_insert_size_distribution",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 4,
    "text_length": 20,
    "sequence_bases": 700,
    "lengths": [100, 100, 200, 300],
    "sha256": "1cc3ddba72c516c92c1b4e49a8cce6e66fe4f916a8034bfc62bccecfb408d5e6"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "source": "numbers",
    "count": 4,
    "excluded_above_max": 0,
    "min": 100,
    "max": 300,
    "mean": 175.0,
    "median": 150.0,
    "stdev": 95.742711,
    "median_absolute_deviation": 50.0,
    "mode": 100,
    "percentiles": {
      "p5": 100.0,
      "p25": 100.0,
      "p75": 225.0,
      "p95": 285.0
    },
    "bin_width": 100,
    "histogram": [
      {
        "start": 100,
        "end": 199,
        "count": 2,
        "fraction": 0.5,
        "cumulative_fraction": 0.5
      },
      {
        "start": 200,
        "end": 299,
        "count": 1,
        "fraction": 0.25,
        "cumulative_fraction": 0.75
      },
      {
        "start": 300,
        "end": 399,
        "count": 1,
        "fraction": 0.25,
        "cumulative_fraction": 1.0
      }
    ]
  }
}
```

## 22) calculate_insert_size_tenx

- **Purpose:** Calculates the fragment size distribution of 10x data from a fragments.tsv table, SAM text or a list of sizes.

- **Method:** Same computation as `calculate_insert_size`; for 10x scATAC `fragments.tsv` (`chrom`, `start`, `end`, `barcode`, `count`) the size is end - start weighted by count. 10x scRNA read 1 holds the barcode and UMI, so scRNA libraries have no read-pair insert size.

- **API:** `dnakit.sequencing_qc.calculate_insert_size_tenx(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input` is a 10x scATAC `fragments.tsv` text or path, SAM text, a list of numbers or a read pair; options as `calculate_insert_size` (`min_overlap`, `max_mismatch_rate`, `max_pairs`, `proper_pairs_only`, `max_insert_size`, `phred_offset`, `max_reads`).

- **Output:** Returns a JSON-serializable envelope with `status` (`"ok"`), `module`, `operation`, `provenance` (`implementation: "native"`, `backend: "local"`), `input_summary`, `warnings` and `result`. `result` holds `source` (`read-pair overlap`, `SAM TLEN`, `fragments (end - start)` or `numbers`), `count`, `excluded_above_max`, `min`, `max`, `mean`, `median`, `stdev`, `median_absolute_deviation`, `mode`, `percentiles` (`p5`, `p25`, `p75`, `p95`) and source-specific counts (`pairs`, `pairs_without_overlap`, `non_proper_pairs_skipped`, ...). Invalid input or options raise `InputFormatError` or `ConfigurationError`.

- **Example code:**

```python
from dnakit.sequencing_qc import calculate_insert_size_tenx

fragments = "chr1\t100\t250\tAAAC\t2\nchr1\t100\t300\tAAAC\t1\n"
result = calculate_insert_size_tenx(fragments)
print(result)
```

- **Example output:**

```json
{
  "module": "sequencing_qc",
  "operation": "calculate_insert_size_tenx",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 2,
    "text_length": 40,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "5a53a6078f92a73b442f3887fee115319bb05572c40c635d00d176af98c64a00"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "source": "fragments (end - start)",
    "count": 3,
    "excluded_above_max": 0,
    "min": 150,
    "max": 200,
    "mean": 166.666667,
    "median": 150.0,
    "stdev": 28.867513,
    "median_absolute_deviation": 0.0,
    "mode": 150,
    "percentiles": {
      "p5": 150.0,
      "p25": 150.0,
      "p75": 175.0,
      "p95": 195.0
    }
  }
}
```

## 23) calculate_phasing_qc

- **Purpose:** Intended to calculate Illumina phasing and pre-phasing percentages from the sequencer's InterOp run metrics.

- **Method:** Not implemented. Phasing and pre-phasing are instrument-side metrics stored in the sequencer's InterOp binaries and cannot be derived from FASTQ reads. The function always returns `status: "not_implemented"` and performs no computation.

- **API:** `dnakit.sequencing_qc.calculate_phasing_qc(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input` is accepted but not used for any computation (a path or any object); `options` is optional and ignored.

- **Output:** Always returns an envelope with `status: "not_implemented"`, `module`, `operation`, `provenance` (`implementation: "not_implemented"`), `input_summary`, `warnings`, `required_dependency` ("Illumina InterOp run-folder files (phasing/prephasing metrics)") and `message` (the reason). There is no `result` key.

- **Example code:**

```python
from dnakit.sequencing_qc import calculate_phasing_qc

result = calculate_phasing_qc("run_folder/InterOp")
print(result)
```

- **Example output:**

```json
{
  "module": "sequencing_qc",
  "operation": "calculate_phasing_qc",
  "provenance": {
    "implementation": "not_implemented",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 1,
    "text_length": 18,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "6dff2776f630ad4ba95c41e73226a139851dee3b5bd6678bac6016b7c1ca929c"
  },
  "warnings": [],
  "status": "not_implemented",
  "required_dependency": "Illumina InterOp run-folder files (phasing/prephasing metrics)",
  "message": "Phasing and pre-phasing are instrument-side metrics stored in the sequencer's InterOp binaries and cannot be derived from FASTQ reads."
}
```

## 24) calculate_phasing_qc_tenx

- **Purpose:** Intended to calculate Illumina phasing and pre-phasing percentages of a 10x library run from the sequencer's InterOp run metrics.

- **Method:** Not implemented. Phasing and pre-phasing are instrument-side metrics stored in the sequencer's InterOp binaries and cannot be derived from FASTQ reads. The function always returns `status: "not_implemented"` and performs no computation.

- **API:** `dnakit.sequencing_qc.calculate_phasing_qc_tenx(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input` is accepted but not used for any computation (a path or any object); `options` is optional and ignored.

- **Output:** Always returns an envelope with `status: "not_implemented"`, `module`, `operation`, `provenance` (`implementation: "not_implemented"`), `input_summary`, `warnings`, `required_dependency` ("Illumina InterOp run-folder files (phasing/prephasing metrics)") and `message` (the reason). There is no `result` key.

- **Example code:**

```python
from dnakit.sequencing_qc import calculate_phasing_qc_tenx

result = calculate_phasing_qc_tenx("run_folder/InterOp")
print(result)
```

- **Example output:**

```json
{
  "module": "sequencing_qc",
  "operation": "calculate_phasing_qc_tenx",
  "provenance": {
    "implementation": "not_implemented",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 1,
    "text_length": 18,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "6dff2776f630ad4ba95c41e73226a139851dee3b5bd6678bac6016b7c1ca929c"
  },
  "warnings": [],
  "status": "not_implemented",
  "required_dependency": "Illumina InterOp run-folder files (phasing/prephasing metrics)",
  "message": "Phasing and pre-phasing are instrument-side metrics stored in the sequencer's InterOp binaries and cannot be derived from FASTQ reads."
}
```

## 25) calculate_qscore_summary

- **Purpose:** Calculates a global Q-score summary: mean and median Q, fractions of bases at Q10-Q40 and read-level Q thresholds.

- **Method:** Base-level statistics are the mean, median, minimum and maximum Phred Q, the Q of the mean error probability and the fraction of bases at or above Q10, Q20, Q30 and Q40. The Q of a read is, by default (`probability`), -10 log10 of its mean error probability (as MinKNOW and NanoPlot report it), or the arithmetic mean of its base Q (`arithmetic`); the number of reads at or above Q7, Q10, Q15, Q20 and Q30 is counted. Pure Python.

- **API:** `dnakit.sequencing_qc.calculate_qscore_summary(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input` is FASTQ text or a file path (`.gz` allowed), a list of reads ((id, sequence, quality) tuples, `{sequence, quality}` mappings, FASTQ record objects or bare sequences), or a pair `{"read1": ..., "read2": ...}` (analysed per mate). Phred qualities are required. Options: `read_q_method` (`probability` default or `arithmetic`), `phred_offset` (33 or 64, default 33), `max_reads` (default 1000000).

- **Output:** Returns a JSON-serializable envelope with `status` (`"ok"`), `module`, `operation`, `provenance` (`implementation: "native"`, `backend: "local"`), `input_summary`, `warnings` and `result`. `result` holds `bases`, `mean_q`, `q_of_mean_error_probability`, `median_q`, `min_q`, `max_q`, `fraction_bases_at_least` (`Q10`, `Q20`, `Q30`, `Q40`), `read_q_method`, `reads`, `mean_read_q`, `median_read_q` and `reads_with_mean_q_at_least` (`Q7`, `Q10`, `Q15`, `Q20`, `Q30`). For a pair the result has `paired: true` and holds these fields separately under `read1` and `read2`. Invalid input or options raise `InputFormatError` or `ConfigurationError`.

- **Example code:**

```python
from dnakit.sequencing_qc import calculate_qscore_summary

def fastq(*records):
    return "".join(f"@{name}\n{seq}\n+\n{qual}\n" for name, seq, qual in records)

reads = fastq(
    ("r1", "ACGTACGTAC", "I" * 10),
    ("r2", "ACGTACGTAC", "I" * 10),
    ("r3", "TTTTTGGGGG", "5" * 10),
    ("r4", "NNACGTACGT", "!!" + "I" * 8),
)
result = calculate_qscore_summary(reads)
print(result)
```

- **Example output:**

```json
{
  "module": "sequencing_qc",
  "operation": "calculate_qscore_summary",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 4,
    "text_length": 112,
    "sequence_bases": 40,
    "lengths": [10, 10, 10, 10],
    "sha256": "74e1d5397e820db4ae7e5d14647765715a6c1e7a262dff425226e4bce9a7169a"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "bases": 40,
    "mean_q": 33.0,
    "q_of_mean_error_probability": 12.79262,
    "median_q": 40.0,
    "min_q": 0,
    "max_q": 40,
    "fraction_bases_at_least": {
      "Q10": 0.95,
      "Q20": 0.95,
      "Q30": 0.7,
      "Q40": 0.7
    },
    "read_q_method": "probability",
    "reads": 4,
    "mean_read_q": 26.746991,
    "median_read_q": 30.0,
    "reads_with_mean_q_at_least": {
      "Q7": 3,
      "Q10": 3,
      "Q15": 3,
      "Q20": 3,
      "Q30": 2
    }
  }
}
```

## 26) calculate_quality

- **Purpose:** Calculates per-position base quality statistics (mean, median, quartiles, 10th/90th percentile) and the per-read mean quality distribution.

- **Method:** For every read position the Phred Q values of the reads covering it give the mean, median, lower and upper quartile, 10th and 90th percentile, minimum and maximum. Percentiles use linear interpolation between order statistics and are rounded to 6 decimals, so floating-point noise such as 6.000000000000001 is not shown. The per-read mean quality is the arithmetic mean of Q (FastQC), summarised as mean, median, mode and a histogram of integer values. Pure Python.

- **API:** `dnakit.sequencing_qc.calculate_quality(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input` is FASTQ text or a file path (`.gz` allowed), a list of reads ((id, sequence, quality) tuples, `{sequence, quality}` mappings, FASTQ record objects or bare sequences), or a pair `{"read1": ..., "read2": ...}` (analysed per mate). Phred qualities are required. Options: `phred_offset` (33 or 64, default 33), `max_reads` (default 1000000).

- **Output:** Returns a JSON-serializable envelope with `status` (`"ok"`), `module`, `operation`, `provenance` (`implementation: "native"`, `backend: "local"`), `input_summary`, `warnings` and `result`. `result` holds `reads_with_quality`, `per_position` (per position `position`, `bases`, `mean`, `median`, `q1`, `q3`, `p10`, `p90`, `min`, `max`) and `per_read_mean_quality` (`mean`, `median`, `mode`, `distribution`). For a pair the result has `paired: true` and holds these fields separately under `read1` and `read2`. Invalid input or options raise `InputFormatError` or `ConfigurationError`.

- **Example code:**

```python
from dnakit.sequencing_qc import calculate_quality

def fastq(*records):
    return "".join(f"@{name}\n{seq}\n+\n{qual}\n" for name, seq, qual in records)

reads = fastq(
    ("r1", "ACGTACGTAC", "I" * 10),
    ("r2", "ACGTACGTAC", "I" * 10),
    ("r3", "TTTTTGGGGG", "5" * 10),
    ("r4", "NNACGTACGT", "!!" + "I" * 8),
)
result = calculate_quality(reads)
print(result)
```

- **Example output:**

```json
{
  "module": "sequencing_qc",
  "operation": "calculate_quality",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 4,
    "text_length": 112,
    "sequence_bases": 40,
    "lengths": [10, 10, 10, 10],
    "sha256": "74e1d5397e820db4ae7e5d14647765715a6c1e7a262dff425226e4bce9a7169a"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "reads_with_quality": 4,
    "per_position": [
      {
        "position": 1,
        "bases": 4,
        "mean": 25.0,
        "median": 30.0,
        "q1": 15.0,
        "q3": 40.0,
        "p10": 6.0,
        "p90": 40.0,
        "min": 0,
        "max": 40
      },
      {
        "position": 2,
        "bases": 4,
        "mean": 25.0,
        "median": 30.0,
        "q1": 15.0,
        "q3": 40.0,
        "p10": 6.0,
        "p90": 40.0,
        "min": 0,
        "max": 40
      },
      "..."
    ],
    "per_read_mean_quality": {
      "mean": 33.0,
      "median": 36.0,
      "mode": 40,
      "distribution": {
        "20": 1,
        "32": 1,
        "40": 2
      }
    }
  }
}
```

## 27) calculate_quality_tenx

- **Purpose:** Calculates per-position quality statistics for the barcode, UMI and RNA regions of 10x reads.

- **Method:** Read 1 is split into barcode and UMI regions and read 2 is the RNA read; each region gets the result of `calculate_quality` (per-position mean, median, quartiles, 10th/90th percentile and the per-read mean quality distribution). The 10x chemistry is not auto-detected: `chemistry` (or the explicit lengths) decides how read 1 is split into cell barcode and UMI.

- **API:** `dnakit.sequencing_qc.calculate_quality_tenx(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input` is a pair `{"read1": <barcode+UMI reads>, "read2": <cDNA reads>}` (FASTQ text, path or read list); a single read set is rejected and read 1 must be at least `barcode_length` + `umi_length` long. Options: `chemistry` (`v1` 14+10, `v2` 16+10, `v3` 16+12 bp barcode+UMI; default `v3`), `barcode_length`, `umi_length`, `phred_offset` (33), `max_reads` (1000000). Phred qualities are required.

- **Output:** Returns a JSON-serializable envelope with `status` (`"ok"`), `module`, `operation`, `provenance` (`implementation: "native"`, `backend: "local"`), `input_summary`, `warnings` and `result`. `result` holds `chemistry`, `barcode_length`, `umi_length`, `read_pairs` and `regions` with `barcode`, `umi` and `rna_read`, each holding `reads_with_quality`, `per_position` and `per_read_mean_quality`. Invalid input or options raise `InputFormatError` or `ConfigurationError`.

- **Example code:**

```python
from dnakit.sequencing_qc import calculate_quality_tenx

def fastq(*records):
    return "".join(f"@{name}\n{seq}\n+\n{qual}\n" for name, seq, qual in records)

barcode = "ACGTACGTACGTACGT"
tenx = {
    "read1": fastq(
        ("p1", barcode + "A" * 12, "I" * 28),
        ("p2", barcode + "A" * 12, "I" * 28),
        ("p3", "T" * 16 + "C" * 12, "I" * 16 + "5" * 12),
    ),
    "read2": fastq(
        ("p1", "GATTACAGATCGGAAGAGCGGTTCAGCAGG", "I" * 30),
        ("p2", "GATTACAGATCGGAAGAGCGGTTCAGCAGG", "I" * 30),
        ("p3", "C" * 30, "5" * 30),
    ),
}
result = calculate_quality_tenx(tenx)
print(result)
```

- **Example output:**

```json
{
  "module": "sequencing_qc",
  "operation": "calculate_quality_tenx",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 1,
    "text_length": 446,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "c2bd8879256656753e784d4d4fe4f7d2de6add20603c60c539095f494a8d364c"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "chemistry": "v3",
    "barcode_length": 16,
    "umi_length": 12,
    "read_pairs": 3,
    "regions": {
      "barcode": {
        "reads_with_quality": 3,
        "per_position": [
          {
            "...": "10 keys"
          },
          {
            "...": "10 keys"
          },
          "..."
        ],
        "per_read_mean_quality": {
          "mean": 40.0,
          "median": 40.0,
          "mode": 40,
          "distribution": {
            "...": "1 keys"
          }
        }
      },
      "umi": {
        "reads_with_quality": 3,
        "per_position": [
          {
            "...": "10 keys"
          },
          {
            "...": "10 keys"
          },
          "..."
        ],
        "per_read_mean_quality": {
          "mean": 33.333333,
          "median": 40.0,
          "mode": 40,
          "distribution": {
            "...": "2 keys"
          }
        }
      },
      "rna_read": {
        "reads_with_quality": 3,
        "per_position": [
          {
            "...": "10 keys"
          },
          {
            "...": "10 keys"
          },
          "..."
        ],
        "per_read_mean_quality": {
          "mean": 33.333333,
          "median": 40.0,
          "mode": 40,
          "distribution": {
            "...": "2 keys"
          }
        }
      }
    }
  }
}
```

## 28) calculate_read_length

- **Purpose:** Calculates the read length distribution: min, max, mean, median, standard deviation, N50, N90, L50 and a histogram.

- **Method:** N50 is the length L such that reads of length at least L contain at least half of all bases (N90 likewise for 90%), and L50 is the number of reads, taken longest first, needed to reach half of the bases. Lengths are histogrammed in bins of `bin_width` (default 1, automatic for wide ranges). Pure Python.

- **API:** `dnakit.sequencing_qc.calculate_read_length(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input` is FASTQ text or a file path (`.gz` allowed), a list of reads ((id, sequence, quality) tuples, `{sequence, quality}` mappings, FASTQ record objects or bare sequences), or a pair `{"read1": ..., "read2": ...}` (analysed per mate). Options: `bin_width` (default 1, automatic for wide ranges), `phred_offset` (33 or 64, default 33), `max_reads` (default 1000000).

- **Output:** Returns a JSON-serializable envelope with `status` (`"ok"`), `module`, `operation`, `provenance` (`implementation: "native"`, `backend: "local"`), `input_summary`, `warnings` and `result`. `result` holds `reads`, `bases`, `min`, `max`, `mean`, `median`, `stdev`, `n50`, `n90`, `l50`, `bin_width` and `histogram` (`start`, `end`, `count`). For a pair the result has `paired: true` and holds these fields separately under `read1` and `read2`. Invalid input or options raise `InputFormatError` or `ConfigurationError`.

- **Example code:**

```python
from dnakit.sequencing_qc import calculate_read_length

reads = ["A" * n for n in (100, 200, 300, 400, 500)]
result = calculate_read_length(reads, options={"bin_width": 100})
print(result)
```

- **Example output:**

```json
{
  "module": "sequencing_qc",
  "operation": "calculate_read_length",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 5,
    "text_length": 1520,
    "sequence_bases": 1500,
    "lengths": [100, 200, 300, 400, 500],
    "sha256": "0212d55422b398acab34c76a9e36b0a6991a4cdefef16a25c8f9b01bea4e23e9"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "reads": 5,
    "bases": 1500,
    "min": 100,
    "max": 500,
    "mean": 300.0,
    "median": 300.0,
    "stdev": 158.113883,
    "n50": 400,
    "n90": 200,
    "l50": 2,
    "bin_width": 100,
    "histogram": [
      {
        "start": 100,
        "end": 199,
        "count": 1
      },
      {
        "start": 200,
        "end": 299,
        "count": 1
      },
      "..."
    ]
  }
}
```

## 29) calculate_read_length_tenx

- **Purpose:** Calculates the read length distribution of read 1 and read 2 of 10x data.

- **Method:** Each mate is summarised with the statistics of `calculate_read_length` (min, max, mean, median, standard deviation, N50, N90, L50, histogram).

- **API:** `dnakit.sequencing_qc.calculate_read_length_tenx(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input` is a pair `{"read1": ..., "read2": ...}` (FASTQ text, path or read list). Options: `bin_width`, `phred_offset`, `max_reads`.

- **Output:** Returns a JSON-serializable envelope with `status` (`"ok"`), `module`, `operation`, `provenance` (`implementation: "native"`, `backend: "local"`), `input_summary`, `warnings` and `result`. `result` holds `read_pairs`, `read1` and `read2`, each with the fields of `calculate_read_length`. Invalid input or options raise `InputFormatError` or `ConfigurationError`.

- **Example code:**

```python
from dnakit.sequencing_qc import calculate_read_length_tenx

def fastq(*records):
    return "".join(f"@{name}\n{seq}\n+\n{qual}\n" for name, seq, qual in records)

barcode = "ACGTACGTACGTACGT"
tenx = {
    "read1": fastq(
        ("p1", barcode + "A" * 12, "I" * 28),
        ("p2", barcode + "A" * 12, "I" * 28),
        ("p3", "T" * 16 + "C" * 12, "I" * 16 + "5" * 12),
    ),
    "read2": fastq(
        ("p1", "GATTACAGATCGGAAGAGCGGTTCAGCAGG", "I" * 30),
        ("p2", "GATTACAGATCGGAAGAGCGGTTCAGCAGG", "I" * 30),
        ("p3", "C" * 30, "5" * 30),
    ),
}
result = calculate_read_length_tenx(tenx)
print(result)
```

- **Example output:**

```json
{
  "module": "sequencing_qc",
  "operation": "calculate_read_length_tenx",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 1,
    "text_length": 446,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "c2bd8879256656753e784d4d4fe4f7d2de6add20603c60c539095f494a8d364c"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "read_pairs": 3,
    "read1": {
      "reads": 3,
      "bases": 84,
      "min": 28,
      "max": 28,
      "mean": 28.0,
      "median": 28.0,
      "stdev": 0.0,
      "n50": 28,
      "n90": 28,
      "l50": 2,
      "bin_width": 1,
      "histogram": [
        {
          "start": 28,
          "end": 28,
          "count": 3
        }
      ]
    },
    "read2": {
      "reads": 3,
      "bases": 90,
      "min": 30,
      "max": 30,
      "mean": 30.0,
      "median": 30.0,
      "stdev": 0.0,
      "n50": 30,
      "n90": 30,
      "l50": 2,
      "bin_width": 1,
      "histogram": [
        {
          "start": 30,
          "end": 30,
          "count": 3
        }
      ]
    }
  }
}
```

## 30) calculate_yield

- **Purpose:** Calculates sequencing yield: reads, bases, mean length, N bases, GC fraction and Q20/Q30 bases.

- **Method:** Counts reads and bases, the mean read length, `N` bases, the GC fraction (G+C over A+C+G+T) and, when qualities exist, the bases at Q20 and Q30 or above and their fractions. With `genome_size` the expected coverage is bases divided by the genome size. Pure Python.

- **API:** `dnakit.sequencing_qc.calculate_yield(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input` is FASTQ text or a file path (`.gz` allowed), a list of reads ((id, sequence, quality) tuples, `{sequence, quality}` mappings, FASTQ record objects or bare sequences), or a pair `{"read1": ..., "read2": ...}` (analysed per mate). Options: `genome_size` (optional; adds `expected_coverage`), `phred_offset` (33 or 64, default 33), `max_reads` (default 1000000).

- **Output:** Returns a JSON-serializable envelope with `status` (`"ok"`), `module`, `operation`, `provenance` (`implementation: "native"`, `backend: "local"`), `input_summary`, `warnings` and `result`. `result` holds `reads`, `bases`, `mean_read_length`, `n_bases`, `gc_fraction`, `q20_bases`, `q30_bases`, `q20_fraction`, `q30_fraction` (`null` without qualities) and, with `genome_size`, `genome_size` and `expected_coverage`. For a pair the result has `paired: true`, `read_pairs`, `total_bases` and these fields under `read1` and `read2`. Invalid input or options raise `InputFormatError` or `ConfigurationError`.

- **Example code:**

```python
from dnakit.sequencing_qc import calculate_yield

def fastq(*records):
    return "".join(f"@{name}\n{seq}\n+\n{qual}\n" for name, seq, qual in records)

reads = fastq(
    ("r1", "ACGTACGTAC", "I" * 10),
    ("r2", "ACGTACGTAC", "I" * 10),
    ("r3", "TTTTTGGGGG", "5" * 10),
    ("r4", "NNACGTACGT", "!!" + "I" * 8),
)
result = calculate_yield(reads, options={"genome_size": 20})
print(result)
```

- **Example output:**

```json
{
  "module": "sequencing_qc",
  "operation": "calculate_yield",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 4,
    "text_length": 112,
    "sequence_bases": 40,
    "lengths": [10, 10, 10, 10],
    "sha256": "74e1d5397e820db4ae7e5d14647765715a6c1e7a262dff425226e4bce9a7169a"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "reads": 4,
    "bases": 40,
    "mean_read_length": 10.0,
    "n_bases": 2,
    "gc_fraction": 0.5,
    "q20_bases": 38,
    "q30_bases": 28,
    "q20_fraction": 0.95,
    "q30_fraction": 0.7,
    "genome_size": 20.0,
    "expected_coverage": 2.0
  }
}
```

## 31) calculate_yield_tenx

- **Purpose:** Calculates the yield of 10x reads per region (barcode, UMI, RNA read) plus barcode counts.

- **Method:** Read 1 is split into barcode and UMI regions and read 2 is the RNA read; each region gets the yield fields of `calculate_yield`. In addition the number of distinct cell barcodes and the median reads per barcode are given and, with a `whitelist`, the fraction of valid barcodes (an exact whitelist entry or a unique whitelist neighbour within `max_mismatches`). The 10x chemistry is not auto-detected: `chemistry` (or the explicit lengths) decides how read 1 is split into cell barcode and UMI.

- **API:** `dnakit.sequencing_qc.calculate_yield_tenx(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input` is a pair `{"read1": <barcode+UMI reads>, "read2": <cDNA reads>}` (FASTQ text, path or read list); a single read set is rejected and read 1 must be at least `barcode_length` + `umi_length` long. Options: `chemistry` (`v1` 14+10, `v2` 16+10, `v3` 16+12 bp barcode+UMI; default `v3`), `barcode_length`, `umi_length`, `phred_offset` (33), `max_reads` (1000000). Extra options: `whitelist` (valid cell barcodes), `max_mismatches` (default 1).

- **Output:** Returns a JSON-serializable envelope with `status` (`"ok"`), `module`, `operation`, `provenance` (`implementation: "native"`, `backend: "local"`), `input_summary`, `warnings` and `result`. `result` holds `chemistry`, `barcode_length`, `umi_length`, `read_pairs`, `regions` (`barcode`, `umi`, `rna_read`, each with the fields of `calculate_yield`), `total_bases_read1_read2`, `distinct_barcodes`, `median_reads_per_barcode` and, with a whitelist, `valid_barcode_fraction` and `valid_barcode_definition`. Invalid input or options raise `InputFormatError` or `ConfigurationError`.

- **Example code:**

```python
from dnakit.sequencing_qc import calculate_yield_tenx

def fastq(*records):
    return "".join(f"@{name}\n{seq}\n+\n{qual}\n" for name, seq, qual in records)

barcode = "ACGTACGTACGTACGT"
tenx = {
    "read1": fastq(
        ("p1", barcode + "A" * 12, "I" * 28),
        ("p2", barcode + "A" * 12, "I" * 28),
        ("p3", "T" * 16 + "C" * 12, "I" * 16 + "5" * 12),
    ),
    "read2": fastq(
        ("p1", "GATTACAGATCGGAAGAGCGGTTCAGCAGG", "I" * 30),
        ("p2", "GATTACAGATCGGAAGAGCGGTTCAGCAGG", "I" * 30),
        ("p3", "C" * 30, "5" * 30),
    ),
}
result = calculate_yield_tenx(tenx, options={"whitelist": [barcode]})
print(result)
```

- **Example output:**

```json
{
  "module": "sequencing_qc",
  "operation": "calculate_yield_tenx",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 1,
    "text_length": 446,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "c2bd8879256656753e784d4d4fe4f7d2de6add20603c60c539095f494a8d364c"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "chemistry": "v3",
    "barcode_length": 16,
    "umi_length": 12,
    "read_pairs": 3,
    "regions": {
      "barcode": {
        "reads": 3,
        "bases": 48,
        "mean_read_length": 16.0,
        "n_bases": 0,
        "gc_fraction": 0.333333,
        "q20_bases": 48,
        "q30_bases": 48,
        "q20_fraction": 1.0,
        "q30_fraction": 1.0
      },
      "umi": {
        "reads": 3,
        "bases": 36,
        "mean_read_length": 12.0,
        "n_bases": 0,
        "gc_fraction": 0.333333,
        "q20_bases": 36,
        "q30_bases": 24,
        "q20_fraction": 1.0,
        "q30_fraction": 0.666667
      },
      "rna_read": {
        "reads": 3,
        "bases": 90,
        "mean_read_length": 30.0,
        "n_bases": 0,
        "gc_fraction": 0.688889,
        "q20_bases": 90,
        "q30_bases": 60,
        "q20_fraction": 1.0,
        "q30_fraction": 0.666667
      }
    },
    "total_bases_read1_read2": 174,
    "distinct_barcodes": 2,
    "median_reads_per_barcode": 1.5,
    "valid_barcode_fraction": 0.666667,
    "valid_barcode_definition": "Barcode equals a whitelist entry or has a unique whitelist neighbour within 1 mismatches."
  }
}
```

## 32) collapse_umi_reads

- **Purpose:** Collapses reads that share a UMI (and read start) into one consensus read per molecule.

- **Method:** Reads are grouped by (UMI, first `start_length` bases) because no alignment position is available, so the read start only approximates the mapping position; with `group_by` `umi` only the UMI is used. With `correct_errors` the UMIs are first merged by the directional network of `correct_umi_errors`. The consensus is a per-position vote weighted by Phred quality, taken over the positions covered by at least half of the group. The consensus quality is clip(sum of Q of agreeing bases - sum of Q of disagreeing bases, 0, 93). Pure Python.

- **API:** `dnakit.sequencing_qc.collapse_umi_reads(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input` is FASTQ reads whose names end in `<umi_separator><UMI>` (umi_tools convention), for example `read1_ACGTACGT`; in the formats FASTQ text or a file path (`.gz` allowed), a list of reads ((id, sequence, quality) tuples, `{sequence, quality}` mappings, FASTQ record objects or bare sequences). A pair mapping is not accepted (`INVALID_INPUT`); for paired data collapse the read that carries the UMI name. Options: `umi_separator` (default `_`), `group_by` (`umi_and_start` default or `umi`), `start_length` (default 10), `correct_errors` (default `False`), `threshold` (Hamming distance for the correction, default 1), `phred_offset` (33 or 64, default 33), `max_reads` (default 1000000).

- **Output:** Returns a JSON-serializable envelope with `status` (`"ok"`), `module`, `operation`, `provenance` (`implementation: "native"`, `backend: "local"`), `input_summary`, `warnings` and `result`. `result` holds `reads_in`, `molecules_out`, `collapse_ratio`, `group_by`, `start_length`, `error_corrected`, `groups` (per group `umi`, `start`, `reads`, `read_names`, `consensus_sequence`, `consensus_quality`) and a `definition` string. Invalid input or options raise `InputFormatError` or `ConfigurationError`.

- **Example code:**

```python
from dnakit.sequencing_qc import collapse_umi_reads

def fastq(*records):
    return "".join(f"@{name}\n{seq}\n+\n{qual}\n" for name, seq, qual in records)

reads = fastq(
    ("r1_AAAA", "ACGTACGT", "IIIIIIII"),
    ("r2_AAAA", "ACGTACGT", "IIIIIIII"),
    ("r3_AAAA", "ACGAACGT", "III5IIII"),
    ("r4_CCCC", "ACGTACGT", "IIIIIIII"),
)
result = collapse_umi_reads(reads, options={"group_by": "umi"})
print(result)
```

- **Example output:**

```json
{
  "module": "sequencing_qc",
  "operation": "collapse_umi_reads",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 4,
    "text_length": 116,
    "sequence_bases": 32,
    "lengths": [8, 8, 8, 8],
    "sha256": "e7ea522ad8a504ec9f45d37aebd7c1e1e0964c3e20b71854138375e0984a97f9"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "reads_in": 4,
    "molecules_out": 2,
    "collapse_ratio": 2.0,
    "group_by": "umi",
    "start_length": null,
    "error_corrected": false,
    "groups": [
      {
        "umi": "AAAA",
        "start": null,
        "reads": 3,
        "read_names": [
          "r1_AAAA",
          "r2_AAAA",
          "..."
        ],
        "consensus_sequence": "ACGTACGT",
        "consensus_quality": "~~~]~~~~"
      },
      {
        "umi": "CCCC",
        "start": null,
        "reads": 1,
        "read_names": [
          "r4_CCCC"
        ],
        "consensus_sequence": "ACGTACGT",
        "consensus_quality": "IIIIIIII"
      }
    ],
    "definition": "Reads are grouped by (UMI, first start_length bases) because no alignment position is available; the consensus is a per-position vote weighted by Phred quality ..."
  }
}
```

## 33) correct_long_reads

- **Purpose:** Intended to correct sequencing errors in long reads by overlap-based consensus (canu/ratatosk style).

- **Method:** Not implemented. Long-read self-correction needs all-vs-all overlap detection and consensus calling that DNAKit does not implement natively. The function always returns `status: "not_implemented"` and performs no computation.

- **API:** `dnakit.sequencing_qc.correct_long_reads(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input` is accepted but not used for any computation (a path or any object); `options` is optional and ignored.

- **Output:** Always returns an envelope with `status: "not_implemented"`, `module`, `operation`, `provenance` (`implementation: "not_implemented"`), `input_summary`, `warnings`, `required_dependency` ("long-read error-correction program (canu, ratatosk or similar)") and `message` (the reason). There is no `result` key.

- **Example code:**

```python
from dnakit.sequencing_qc import correct_long_reads

result = correct_long_reads("long_reads.fastq")
print(result)
```

- **Example output:**

```json
{
  "module": "sequencing_qc",
  "operation": "correct_long_reads",
  "provenance": {
    "implementation": "not_implemented",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 1,
    "text_length": 16,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "b8567bebd311c7484fd0d8dbde61c467727df0b1d64337a7ef96d0dac8922bf8"
  },
  "warnings": [],
  "status": "not_implemented",
  "required_dependency": "long-read error-correction program (canu, ratatosk or similar)",
  "message": "Long-read self-correction needs all-vs-all overlap detection and consensus calling that DNAKit does not implement natively."
}
```

## 34) correct_umi_errors

- **Purpose:** Corrects UMI sequencing errors by clustering UMIs on a Hamming-distance network (UMI-tools methods).

- **Method:** UMIs within Hamming distance `threshold` are linked and grouped by one of three UMI-tools methods. `directional` (default): an edge a -> b exists when count(a) >= 2 x count(b) - 1, and groups are grown by breadth-first search from the most abundant UMI. `adjacency`: the most abundant UMI plus its neighbours, repeated within each connected component. `cluster`: connected components. The most abundant UMI of each group is its representative and the other members are corrected to it. Pure Python.

- **API:** `dnakit.sequencing_qc.correct_umi_errors(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input` is a list of UMI strings, a mapping UMI -> read count, or reads whose names end in `_<UMI>` (at most 5000 distinct UMIs). Options: `method` (`directional` default, `adjacency` or `cluster`), `threshold` (Hamming distance, default 1), `umi_separator` (default `_`), `phred_offset` (33 or 64, default 33), `max_reads` (default 1000000).

- **Output:** Returns a JSON-serializable envelope with `status` (`"ok"`), `module`, `operation`, `provenance` (`implementation: "native"`, `backend: "local"`), `input_summary`, `warnings` and `result`. `result` holds `method`, `threshold`, `distinct_umis_before`, `distinct_umis_after`, `reads`, `groups` (per group `representative`, `members`, `reads`, `member_counts`), `corrections` (corrected UMI -> representative) and a `definition` string. Invalid input or options raise `InputFormatError` or `ConfigurationError`.

- **Example code:**

```python
from dnakit.sequencing_qc import correct_umi_errors

counts = {"AAAAAA": 100, "AAAAAT": 1, "CCCCCC": 50, "CCCCCG": 30}
result = correct_umi_errors(counts)
print(result)
```

- **Example output:**

```json
{
  "module": "sequencing_qc",
  "operation": "correct_umi_errors",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 1,
    "text_length": 56,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "3fefe1abde49ead8296fe12856ea08439dda7deb9e9f8853bab0cb087b034aba"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "method": "directional",
    "threshold": 1,
    "distinct_umis_before": 4,
    "distinct_umis_after": 3,
    "reads": 181,
    "groups": [
      {
        "representative": "AAAAAA",
        "members": [
          "AAAAAT"
        ],
        "reads": 101,
        "member_counts": {
          "AAAAAA": 100,
          "AAAAAT": 1
        }
      },
      {
        "representative": "CCCCCC",
        "members": [],
        "reads": 50,
        "member_counts": {
          "CCCCCC": 50
        }
      },
      "..."
    ],
    "corrections": {
      "AAAAAT": "AAAAAA"
    },
    "definition": "UMI-tools networks over Hamming distance <= threshold. directional: edge a->b when count(a) >= 2*count(b)-1, groups grown by breadth-first search from the most ..."
  }
}
```

## 35) demultiplex_barcoded_reads

- **Purpose:** Demultiplexes reads into per-sample read sets by barcode, trimming the barcode when it is part of the read.

- **Method:** Barcode assignment is that of `assign_sample_barcodes` (Hamming distance with `max_mismatches` and `min_distance_gap`). Assigned reads are returned per sample; for `read_start` and `read_end` barcodes the barcode bases (and their qualities) are trimmed when `trim_barcode` is true. Reads without an unambiguous match go to `unassigned` and are not trimmed. With `output_directory` one `<sample>.fastq` per sample is written, and an existing file is never overwritten. Pure Python.

- **API:** `dnakit.sequencing_qc.demultiplex_barcoded_reads(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input` and options as `assign_sample_barcodes` (`samples` required, `barcode_location`, `barcode_length`, `max_mismatches`, `min_distance_gap`; a pair mapping is not accepted and raises `INVALID_INPUT`; for paired data demultiplex the read that carries the barcode), plus `trim_barcode` (default `True` for `read_start` and `read_end`), `output_directory` (optional existing directory) and `max_listed` (reads returned per sample, default 1000).

- **Output:** Returns a JSON-serializable envelope with `status` (`"ok"`), `module`, `operation`, `provenance` (`implementation: "native"`, `backend: "local"`), `input_summary`, `warnings` and `result`. `result` holds the assignment summary of `assign_sample_barcodes` (`reads`, `reads_per_sample`, `assigned`, `unassigned`, `ambiguous`, `assigned_fraction`, ...), `barcode_trimmed`, `samples` (per sample and `unassigned`: `reads` and `records` with `id`, `sequence`, `quality`) and `written_files` (`null` unless `output_directory` was given). Invalid input or options raise `InputFormatError` or `ConfigurationError`.

- **Example code:**

```python
from dnakit.sequencing_qc import demultiplex_barcoded_reads

def fastq(*records):
    return "".join(f"@{name}\n{seq}\n+\n{qual}\n" for name, seq, qual in records)

reads = fastq(
    ("a", "ACGTAAAA", "IIIIJJJJ"),
    ("b", "TTTTCCCC", "IIIIJJJJ"),
    ("c", "GGGGCCCC", "IIIIJJJJ"),
)
result = demultiplex_barcoded_reads(reads, options={"samples": {"S1": "ACGT", "S2": "TTTT"}, "barcode_length": 4})
print(result)
```

- **Example output:**

```json
{
  "module": "sequencing_qc",
  "operation": "demultiplex_barcoded_reads",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 3,
    "text_length": 69,
    "sequence_bases": 24,
    "lengths": [8, 8, 8],
    "sha256": "554e09db1fbcb9792fc57c3561da72f69cb5cf168c635f51759ab0ce0ac285a1"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "reads": 3,
    "barcode_location": "read_start",
    "max_mismatches": 1,
    "min_distance_gap": 1,
    "reads_per_sample": {
      "S1": 1,
      "S2": 1
    },
    "assigned": 2,
    "unassigned": 1,
    "ambiguous": 0,
    "assigned_fraction": 0.666667,
    "min_distance_between_sample_barcodes": 3,
    "barcode_trimmed": true,
    "samples": {
      "S1": {
        "reads": 1,
        "records": [
          {
            "...": "3 keys"
          }
        ]
      },
      "S2": {
        "reads": 1,
        "records": [
          {
            "...": "3 keys"
          }
        ]
      },
      "unassigned": {
        "reads": 1,
        "records": [
          {
            "...": "3 keys"
          }
        ]
      }
    },
    "written_files": null
  }
}
```

## 36) detect_cross_sample_contamination

- **Purpose:** Flags sample pairs that share an unexpectedly large part of their k-mers (possible cross-sample contamination).

- **Method:** For every sample the set of canonical k-mers is built, and containment[a][b] = |K(a) intersect K(b)| / |K(a)| is computed for each ordered pair. Pairs with containment at or above `threshold` are flagged. Heuristic: it is only meaningful when the samples are expected to be unrelated (different organisms or individuals); related samples share k-mers legitimately. Pure Python.

- **API:** `dnakit.sequencing_qc.detect_cross_sample_contamination(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input` is a mapping sample name -> reads (at least two samples; each in the formats FASTQ text or a file path (`.gz` allowed), a list of reads ((id, sequence, quality) tuples, `{sequence, quality}` mappings, FASTQ record objects or bare sequences)). Options: `kmer_size` (default 21), `threshold` (containment in [0, 1], default 0.1), `phred_offset` (33 or 64, default 33), `max_reads` (default 1000000).

- **Output:** Returns a JSON-serializable envelope with `status` (`"ok"`), `module`, `operation`, `provenance` (`implementation: "native"`, `backend: "local"`), `input_summary`, `warnings` and `result`. `result` holds `samples` (sample -> number of distinct k-mers), `kmer_size`, `threshold`, `containment` (nested mapping a -> b -> fraction), `flagged_pairs` (`sample`, `found_in`, `containment`, sorted by decreasing containment) and a `definition` string. Invalid input or options raise `InputFormatError` or `ConfigurationError`.

- **Example code:**

```python
from dnakit.sequencing_qc import detect_cross_sample_contamination

host = "ACGTTGCAAGGCTTAACCGGATCGATTCGAGGCTATTGCCAGTCAAGGTACCTGATGCCGTAAGCTTGGAACCTAGGTCCA"
other = "TTGACCGGTAACCTTGGAACCGTTAGGCCATTGACCAAGGTTCCAAGGTTCCGGAATTCCGGAACCTTGGAAGGCC"
samples = {"A": [host[:40]], "B": [host[20:60]], "C": [other[:40]]}
result = detect_cross_sample_contamination(samples, options={"kmer_size": 12, "threshold": 0.2})
print(result)
```

- **Example output:**

```json
{
  "module": "sequencing_qc",
  "operation": "detect_cross_sample_contamination",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 1,
    "text_length": 153,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "8435281fcc502972cb1960fbe982a9f256340d91c17229f1e34dadc196a800c1"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "samples": {
      "A": 29,
      "B": 29,
      "C": 29
    },
    "kmer_size": 12,
    "threshold": 0.2,
    "containment": {
      "A": {
        "B": 0.310345,
        "C": 0.0
      },
      "B": {
        "A": 0.310345,
        "C": 0.0
      },
      "C": {
        "A": 0.0,
        "B": 0.0
      }
    },
    "flagged_pairs": [
      {
        "sample": "A",
        "found_in": "B",
        "containment": 0.310345
      },
      {
        "sample": "B",
        "found_in": "A",
        "containment": 0.310345
      }
    ],
    "definition": "containment[a][b] = |K(a) intersect K(b)| / |K(a)| over canonical k-mers. Heuristic: samples that are expected to be biologically unrelated should share almost ..."
  }
}
```

## 37) detect_platform

- **Purpose:** Detects the sequencing platform from the read-name format (rule-based).

- **Method:** The name of each read is matched against known layouts: Illumina (instrument:run:flowcell:lane:tile:x:y), Oxford Nanopore (`runid=`, `ch=`, UUID), PacBio (movie/zmw/ccs), MGI/BGI (L#C#R#) and Ion Torrent; SRA-renamed reads carry no platform information. The platform is the majority layout when it covers at least half of the reads, otherwise `unknown`. The instrument model is not inferred. Pure Python.

- **API:** `dnakit.sequencing_qc.detect_platform(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input` is FASTQ text or a file path (`.gz` allowed), a list of reads ((id, sequence, quality) tuples, `{sequence, quality}` mappings, FASTQ record objects or bare sequences), or a pair `{"read1": ..., "read2": ...}` (analysed per mate). Options: `phred_offset` (33 or 64, default 33), `max_reads` (default 1000000).

- **Output:** Returns a JSON-serializable envelope with `status` (`"ok"`), `module`, `operation`, `provenance` (`implementation: "native"`, `backend: "local"`), `input_summary`, `warnings` and `result`. `result` holds `platform` (`illumina`, `oxford_nanopore`, `pacbio`, `mgi_dnbseq`, `ion_torrent`, `sra_renamed` or `unknown`), `confidence_fraction`, `read_name_votes`, `read_length` (`min`, `median`, `max`, `uniform`), `evidence` (for Illumina `instrument_id`, `run_number`, `flowcell`, `lane`) and a `definition` string. For a pair the result has `paired: true` and holds these fields separately under `read1` and `read2`. Invalid input or options raise `InputFormatError` or `ConfigurationError`.

- **Example code:**

```python
from dnakit.sequencing_qc import detect_platform

reads = "@A00123:45:HXYZ123DSXX:1:1101:1000:2000 1:N:0:ACGT\nACGT\n+\nIIII\n"
result = detect_platform(reads)
print(result)
```

- **Example output:**

```json
{
  "module": "sequencing_qc",
  "operation": "detect_platform",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 1,
    "text_length": 63,
    "sequence_bases": 4,
    "lengths": [4],
    "sha256": "104d5bfb8cfd56a46818297271296941ac13efeb1f9873b9e9ad65abc7e35c8e"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "platform": "illumina",
    "confidence_fraction": 1.0,
    "read_name_votes": {
      "illumina": 1
    },
    "read_length": {
      "min": 4,
      "median": 4.0,
      "max": 4,
      "uniform": true
    },
    "evidence": {
      "instrument_id": "A00123",
      "run_number": "45",
      "flowcell": "HXYZ123DSXX",
      "lane": "1"
    },
    "definition": "Rule-based: the read-name (header) format of each read is matched against known layouts (Illumina instrument:run:flowcell:lane:tile:x:y, Oxford Nanopore runid=/..."
  }
}
```

## 38) detect_platform_tenx

- **Purpose:** Detects the sequencing platform of 10x reads from read-2 names and infers the Chromium barcode+UMI layout from the read-1 length.

- **Method:** The platform is detected from the read-2 names as in `detect_platform`. Only the modal read-1 length is used (the read-1 sequences are not compared with a chemistry or whitelist): 24, 26 or 28 maps to a 14+10, 16+10 or 16+12 bp barcode+UMI layout and the compatible Chromium chemistries. This is a length heuristic: Cell Ranger's whitelist-based chemistry autodetection is not performed, and other read-1 lengths give no layout.

- **API:** `dnakit.sequencing_qc.detect_platform_tenx(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input` is a pair `{"read1": <barcode+UMI reads>, "read2": <cDNA reads>}` (FASTQ text, path or read list); a single read set is rejected. Options: `phred_offset` (33), `max_reads` (1000000). The chemistry options of the other 10x functions are not used here: the layout is inferred from the read-1 length.

- **Output:** Returns a JSON-serializable envelope with `status` (`"ok"`), `module`, `operation`, `provenance` (`implementation: "native"`, `backend: "local"`), `input_summary`, `warnings` and `result`. `result` holds the fields of `detect_platform` for read 2 (`platform`, `confidence_fraction`, `read_name_votes`, `read_length`, `evidence`, `definition`) and `tenx` (`read1_modal_length`, `layout`, `compatible_chemistries`, `note`; `layout` and `compatible_chemistries` are `null` for an unknown length). Invalid input or options raise `InputFormatError` or `ConfigurationError`.

- **Example code:**

```python
from dnakit.sequencing_qc import detect_platform_tenx

def fastq(*records):
    return "".join(f"@{name}\n{seq}\n+\n{qual}\n" for name, seq, qual in records)

name = "A00123:45:HXYZ123DSXX:1:1101:1000:2000"
tenx = {
    "read1": fastq((name, "ACGTACGTACGTACGT" + "A" * 12, "I" * 28)),
    "read2": fastq((name, "ACGT" * 5, "I" * 20)),
}
result = detect_platform_tenx(tenx)
print(result)
```

- **Example output:**

```json
{
  "module": "sequencing_qc",
  "operation": "detect_platform_tenx",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 1,
    "text_length": 218,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "1b27ae92e3ff96f5c0a9fe0740c51141c29144ec1cc927301189bc3b8dc3d586"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "platform": "illumina",
    "confidence_fraction": 1.0,
    "read_name_votes": {
      "illumina": 1
    },
    "read_length": {
      "min": 20,
      "median": 20.0,
      "max": 20,
      "uniform": true
    },
    "evidence": {
      "instrument_id": "A00123",
      "run_number": "45",
      "flowcell": "HXYZ123DSXX",
      "lane": "1"
    },
    "definition": "Rule-based: the read-name (header) format of each read is matched against known layouts (Illumina instrument:run:flowcell:lane:tile:x:y, Oxford Nanopore runid=/...",
    "tenx": {
      "read1_modal_length": 28,
      "layout": "16 bp barcode + 12 bp UMI",
      "compatible_chemistries": "Chromium 3' v3 / v3.1 / 5' HT",
      "note": "Chemistry is inferred only from the modal read-1 length (24/26/28). Cell Ranger's whitelist-based autodetection needs the 10x barcode whitelists and is not perf..."
    }
  }
}
```

## 39) estimate_library_complexity

- **Purpose:** Estimates library complexity (number of distinct molecules) with the Lander-Waterman / Picard model.

- **Method:** With N total reads and C distinct reads the library size X solves C / X = 1 - exp(-N / X) (Picard EstimateLibraryComplexity), found numerically; the unique reads expected at 0.5x, 1x, 2x, 5x and 10x of the current depth are predicted as X (1 - exp(-k N / X)). X is unbounded (`null`) when every read is distinct. When reads are given, the distinct set is defined by exact-sequence duplicates (pairs use both mates). Pure Python.

- **API:** `dnakit.sequencing_qc.estimate_library_complexity(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input` is reads (FASTQ text or a file path (`.gz` allowed), a list of reads ((id, sequence, quality) tuples, `{sequence, quality}` mappings, FASTQ record objects or bare sequences), or a pair) or a mapping `{"total_reads": N, "unique_reads": C}` (also accepted in `options`). Options: `prefix_length`, `reverse_complement_aware`, `phred_offset` (33 or 64, default 33), `max_reads` (default 1000000).

- **Output:** Returns a JSON-serializable envelope with `status` (`"ok"`), `module`, `operation`, `provenance` (`implementation: "native"`, `backend: "local"`), `input_summary`, `warnings` and `result`. `result` holds `source` (`counts` or reads), `total_reads`, `unique_reads`, `duplicate_rate`, `estimated_library_size` (`null` if unbounded), `predicted_unique_reads` (`0.5x`, `1x`, `2x`, `5x`, `10x`) and a `definition` string. Invalid input or options raise `InputFormatError` or `ConfigurationError`.

- **Example code:**

```python
from dnakit.sequencing_qc import estimate_library_complexity

result = estimate_library_complexity({"total_reads": 1_000_000, "unique_reads": 632_121})
print(result)
```

- **Example output:**

```json
{
  "module": "sequencing_qc",
  "operation": "estimate_library_complexity",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 1,
    "text_length": 48,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "34f4ef63ded2460579d54590e66569cec3c6acf812aa7ca34e25176a55cfe90a"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "source": "counts",
    "total_reads": 1000000,
    "unique_reads": 632121,
    "duplicate_rate": 0.367879,
    "estimated_library_size": 1000001.669581,
    "predicted_unique_reads": {
      "0.5x": 393469.49089,
      "1x": 632121.0,
      "2x": 864665.708484,
      "5x": 993263.655084,
      "10x": 999956.268817
    },
    "definition": "Lander-Waterman / Picard EstimateLibraryComplexity: find X with unique/X = 1 - exp(-total/X)."
  }
}
```

## 40) estimate_read_error_model

- **Purpose:** Estimates a read error model: observed mismatch rates against a reference, or the Phred-derived expected error profile.

- **Method:** With a `reference`, each read is placed on the reference by k-mer seed voting (both strands, the top 3 diagonals) and scored by gap-free Hamming distance; reads above `max_mismatch_fraction` are not counted. The result gives the overall and per-position mismatch rate, substitution counts and, where qualities exist, the empirical error rate per reported Q. Insertions and deletions are not modelled, so indel-rich reads are excluded or inflate the substitution rate. Without a reference the model is the quality-derived profile of `calculate_error_profile` (qualities required). Pure Python.

- **API:** `dnakit.sequencing_qc.estimate_read_error_model(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input` is FASTQ text or a file path (`.gz` allowed), a list of reads ((id, sequence, quality) tuples, `{sequence, quality}` mappings, FASTQ record objects or bare sequences). A pair mapping is not accepted (`INVALID_INPUT`); for paired data run the function once per mate. Options: `reference` (sequence or FASTA text up to 5 Mb; without it the model is the quality-derived profile), `seed_length` (default 11), `max_mismatch_fraction` (default 0.2), `phred_offset` (33 or 64, default 33), `max_reads` (default 1000000).

- **Output:** Returns a JSON-serializable envelope with `status` (`"ok"`), `module`, `operation`, `provenance` (`implementation: "native"`, `backend: "local"`), `input_summary`, `warnings` and `result`. `result` holds `model` with, for a reference, `source` (`ungapped_reference_alignment`), `reads`, `reads_aligned`, `bases_compared`, `mismatches`, `overall_mismatch_rate`, `per_position_mismatch_rate`, `substitution_counts` (for example `A>T`), `empirical_error_by_reported_q` and a `definition`; without a reference the fields of `calculate_error_profile` plus a `note`. Invalid input or options raise `InputFormatError` or `ConfigurationError`.

- **Example code:**

```python
from dnakit.sequencing_qc import estimate_read_error_model

def fastq(*records):
    return "".join(f"@{name}\n{seq}\n+\n{qual}\n" for name, seq, qual in records)

host = "ACGTTGCAAGGCTTAACCGGATCGATTCGAGGCTATTGCCAGTCAAGGTACCTGATGCCGTAAGCTTGGAACCTAGGTCCA"
other = "TTGACCGGTAACCTTGGAACCGTTAGGCCATTGACCAAGGTTCCAAGGTTCCGGAATTCCGGAACCTTGGAAGGCC"
def reverse_complement(seq):
    return seq.translate(str.maketrans("ACGT", "TGCA"))[::-1]

mutate = {"A": "C", "C": "G", "G": "T", "T": "A"}
read1 = host[:10] + mutate[host[10]] + host[11:30]
read2 = host[20:50]
original = reverse_complement(host[40:70])
read3 = mutate[original[0]] + original[1:]
reads = fastq(("r1", read1, "I" * 30), ("r2", read2, "5" * 30), ("r3", read3, "I" * 30))
result = estimate_read_error_model(reads, options={"reference": host})
print(result)
```

- **Example output:**

```json
{
  "module": "sequencing_qc",
  "operation": "estimate_read_error_model",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 3,
    "text_length": 204,
    "sequence_bases": 90,
    "lengths": [30, 30, 30],
    "sha256": "7a26ebdc2713a3ab7bca25d5b393a92b5aa85a9975f53b35f4538570801cd5bc"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "model": {
      "source": "ungapped_reference_alignment",
      "reads": 3,
      "reads_aligned": 3,
      "bases_compared": 90,
      "mismatches": 2,
      "overall_mismatch_rate": 0.022222,
      "per_position_mismatch_rate": [
        {
          "position": 1,
          "bases": 3,
          "mismatch_rate": 0.333333
        },
        {
          "position": 2,
          "bases": 3,
          "mismatch_rate": 0.0
        },
        "..."
      ],
      "substitution_counts": {
        "A>T": 1,
        "G>T": 1
      },
      "empirical_error_by_reported_q": {
        "20": {
          "bases": 30,
          "empirical_error_rate": 0.0,
          "reported_error_rate": 0.01
        },
        "40": {
          "bases": 60,
          "empirical_error_rate": 0.033333,
          "reported_error_rate": 0.0001
        }
      },
      "definition": "Each read is placed on the reference by k-mer seed voting (both strands, top 3 diagonals) and scored by Hamming distance without gaps; reads above max_mismatch_..."
    }
  }
}
```

## 41) extract_umi_sequences

- **Purpose:** Extracts UMI sequences from reads and appends them to the read names (umi_tools convention).

- **Method:** The first (`5prime`) or last (`3prime`) `umi_length` bases are cut from each read and appended to its name as `<name><umi_separator><UMI>`. With `pattern` (5' end only) the positions are given per base: `N` = UMI base, `C` = cell-barcode base, `X` = discarded base, for example `CCCCCCNNNNXX`; the remainder of the read is kept. The trimmed sequence and quality are returned. Reads that are not longer than the extracted region are rejected. Pure Python.

- **API:** `dnakit.sequencing_qc.extract_umi_sequences(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input` is FASTQ text or a file path (`.gz` allowed), a list of reads ((id, sequence, quality) tuples, `{sequence, quality}` mappings, FASTQ record objects or bare sequences). A pair mapping is not accepted (`INVALID_INPUT`); for paired data pass the read that carries the UMI (usually read 1). Options: `umi_length` (default 8) with `umi_location` (`5prime` default or `3prime`), or `pattern` (5' only), `umi_separator` (default `_`), `append_to_name` (default `True`), `phred_offset` (33 or 64, default 33), `max_reads` (default 1000000).

- **Output:** Returns a JSON-serializable envelope with `status` (`"ok"`), `module`, `operation`, `provenance` (`implementation: "native"`, `backend: "local"`), `input_summary`, `warnings` and `result`. `result` holds `reads`, `umi_location`, `umi_region_length`, `distinct_umis`, `records` (per read `id`, `umi`, `cell_barcode`, `sequence`, `quality`) and a `definition` string. Invalid input or options raise `InputFormatError` or `ConfigurationError`.

- **Example code:**

```python
from dnakit.sequencing_qc import extract_umi_sequences

def fastq(*records):
    return "".join(f"@{name}\n{seq}\n+\n{qual}\n" for name, seq, qual in records)

reads = fastq(("r1", "ACGTAAAATTTT", "ABCDEFGHIJKL"))
result = extract_umi_sequences(reads, options={"umi_length": 4})
print(result)
```

- **Example output:**

```json
{
  "module": "sequencing_qc",
  "operation": "extract_umi_sequences",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 1,
    "text_length": 32,
    "sequence_bases": 12,
    "lengths": [12],
    "sha256": "b3a517a420b416cdb6f02aafbee98222c72950c3c2b06697d60edf9b11a7e7fd"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "reads": 1,
    "umi_location": "5prime",
    "umi_region_length": 4,
    "distinct_umis": 1,
    "records": [
      {
        "id": "r1_ACGT",
        "umi": "ACGT",
        "cell_barcode": null,
        "sequence": "AAAATTTT",
        "quality": "EFGHIJKL"
      }
    ],
    "definition": "The first (5prime) or last (3prime) umi_length bases, or the N/C/X pattern positions at the 5' end, are removed from each read and appended to its name as '<nam..."
  }
}
```

## 42) generate_qc_report

- **Purpose:** Generates a FastQC-style QC report: yield, length, quality, GC, adapter, duplication and platform sections with pass/warn/fail flags.

- **Method:** The section `reads` collects the results of `calculate_yield`, `calculate_read_length`, `detect_platform`, `calculate_gc_bias`, `calculate_adapter_content`, `calculate_duplication`, `calculate_qscore_summary` and `calculate_quality` in compact form. Six flags follow the FastQC thresholds listed in the result (`flag_definitions`): per-base quality, per-sequence quality, per-base N content, adapter content, sequence duplication and GC-content shape (`not_available` when it cannot be computed). A plain-text report (`report_text`) is generated and written to a file only when `output_path` is given. Pure Python; it does not call FastQC.

- **API:** `dnakit.sequencing_qc.generate_qc_report(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input` is FASTQ text or a file path (`.gz` allowed), a list of reads ((id, sequence, quality) tuples, `{sequence, quality}` mappings, FASTQ record objects or bare sequences), or a pair `{"read1": ..., "read2": ...}` (analysed per mate). Options: `adapters` (as `calculate_adapter_content`), `output_path` (optional `.txt` or `.md` file; an existing file is not overwritten unless `overwrite=True`), `phred_offset` (33 or 64, default 33), `max_reads` (default 1000000).

- **Output:** Returns a JSON-serializable envelope with `status` (`"ok"`), `module`, `operation`, `provenance` (`implementation: "native"`, `backend: "local"`), `input_summary`, `warnings` and `result`. `result` holds `paired`, `sections` (`reads`, or `read1` and `read2` for a pair, each with `yield`, `read_length`, `platform`, `gc_content`, `adapter_content`, `duplication`, `max_n_percent_at_a_position`, `qscore`, `per_base_quality`, `per_read_mean_quality` and `flags`), `pair_duplication` (pairs only), `flag_definitions`, `report_text` and `output_path` (`null` unless a file was written). Invalid input or options raise `InputFormatError` or `ConfigurationError`.

- **Example code:**

```python
from dnakit.sequencing_qc import generate_qc_report

def fastq(*records):
    return "".join(f"@{name}\n{seq}\n+\n{qual}\n" for name, seq, qual in records)

reads = fastq(
    ("r1", "ACGTACGTAC", "I" * 10),
    ("r2", "ACGTACGTAC", "I" * 10),
    ("r3", "TTTTTGGGGG", "5" * 10),
    ("r4", "NNACGTACGT", "!!" + "I" * 8),
)
result = generate_qc_report(reads)
print(result)
```

- **Example output:**

```json
{
  "module": "sequencing_qc",
  "operation": "generate_qc_report",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 4,
    "text_length": 112,
    "sequence_bases": 40,
    "lengths": [10, 10, 10, 10],
    "sha256": "74e1d5397e820db4ae7e5d14647765715a6c1e7a262dff425226e4bce9a7169a"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "paired": false,
    "sections": {
      "reads": {
        "yield": {
          "reads": 4,
          "bases": 40,
          "mean_read_length": 10.0,
          "n_bases": 2,
          "gc_fraction": 0.5,
          "q20_bases": 38,
          "q30_bases": 28,
          "q20_fraction": 0.95,
          "q30_fraction": 0.7
        },
        "read_length": {
          "reads": 4,
          "bases": 40,
          "min": 10,
          "max": 10,
          "mean": 10.0,
          "median": 10.0,
          "stdev": 0.0,
          "n50": 10,
          "n90": 10,
          "l50": 2,
          "bin_width": 1
        },
        "platform": {
          "platform": "unknown",
          "confidence_fraction": 1.0,
          "read_name_votes": {
            "...": "1 keys"
          },
          "read_length": {
            "...": "4 keys"
          },
          "evidence": {
            "...": "0 keys"
          }
        },
        "gc_content": {
          "mean_gc_percent": 50.0,
          "stdev_gc_percent": 0.0,
          "mode_gc_percent": 50,
          "deviation_from_normal_percent": null
        },
        "adapter_content": {
          "fraction_reads_with_any_adapter": 0.0,
          "max_cumulative_percent": {
            "...": "5 keys"
          }
        },
        "duplication": {
          "distinct_sequences": 3,
          "percent_remaining_after_deduplication": 75.0
        },
        "max_n_percent_at_a_position": 25.0,
        "qscore": {
          "bases": 40,
          "mean_q": 33.0,
          "q_of_mean_error_probability": 12.79262,
          "median_q": 40.0,
          "min_q": 0,
          "max_q": 40,
          "fraction_bases_at_least": {
            "...": "4 keys"
          },
          "read_q_method": "probability",
          "reads": 4,
          "mean_read_q": 26.746991,
          "median_read_q": 30.0,
          "reads_with_mean_q_at_least": {
            "...": "5 keys"
          }
        },
        "per_base_quality": {
          "positions": 10,
          "min_median": 30.0,
          "min_lower_quartile": 15.0
        },
        "per_read_mean_quality": {
          "mean": 33.0,
          "median": 36.0,
          "mode": 40
        },
        "flags": {
          "per_base_quality": "pass",
          "per_sequence_quality": "pass",
          "per_base_n_content": "fail",
          "adapter_content": "pass",
          "sequence_duplication": "warn",
          "gc_content_shape": "not_available"
        }
      }
    },
    "flag_definitions": "FastQC thresholds: per-base quality warn if lower quartile < 10 or median < 25 (fail < 5 / < 20) at any position; per-sequence quality warn if modal read mean Q...",
    "report_text": "DNAKit sequencing QC report\n===========================\n\n[reads]\nreads: 4  bases: 40  mean length: 10.0\nQ20 fraction: 0.95  Q30 fraction: 0.7\nmean GC%: 50.0  pl...",
    "output_path": null
  }
}
```

## 43) generate_qc_report_tenx

- **Purpose:** Generates a QC report for 10x Chromium reads: yield and basecalling per region, barcode/UMI validity, a duplication proxy and FastQC-style flags for the RNA read.

- **Method:** Barcode, UMI and RNA-read regions get the yield and basecalling metrics of `calculate_yield` and `calculate_basecalling_qc`; the RNA read also gets the compact FastQC-style section and flags of `generate_qc_report` (the flags are present when Phred qualities exist). The report adds the barcode/UMI report of `calculate_barcode_quality_tenx`, the duplication proxy of `calculate_duplication_tenx` and the yield summary of `calculate_yield_tenx`. A text report is returned and written to a file only when `output_path` is given. The 10x chemistry is not auto-detected: `chemistry` (or the explicit lengths) decides how read 1 is split into cell barcode and UMI.

- **API:** `dnakit.sequencing_qc.generate_qc_report_tenx(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input` is a pair `{"read1": <barcode+UMI reads>, "read2": <cDNA reads>}` (FASTQ text, path or read list); a single read set is rejected and read 1 must be at least `barcode_length` + `umi_length` long. Options: `chemistry` (`v1` 14+10, `v2` 16+10, `v3` 16+12 bp barcode+UMI; default `v3`), `barcode_length`, `umi_length`, `phred_offset` (33), `max_reads` (1000000). Extra options: `whitelist`, `adapters`, `output_path` (optional `.txt` or `.md`; no overwrite unless `overwrite=True`).

- **Output:** Returns a JSON-serializable envelope with `status` (`"ok"`), `module`, `operation`, `provenance` (`implementation: "native"`, `backend: "local"`), `input_summary`, `warnings` and `result`. `result` holds `chemistry`, `barcode_length`, `umi_length`, `read_pairs`, `sections` (`barcode`, `umi` and `rna_read` with `yield` and `basecalling`; `rna_read` also `compact`), `barcode_and_umi`, `duplication`, `yield`, `report_text` and `output_path` (`null` unless a file was written). Invalid input or options raise `InputFormatError` or `ConfigurationError`.

- **Example code:**

```python
from dnakit.sequencing_qc import generate_qc_report_tenx

def fastq(*records):
    return "".join(f"@{name}\n{seq}\n+\n{qual}\n" for name, seq, qual in records)

barcode = "ACGTACGTACGTACGT"
tenx = {
    "read1": fastq(
        ("p1", barcode + "A" * 12, "I" * 28),
        ("p2", barcode + "A" * 12, "I" * 28),
        ("p3", "T" * 16 + "C" * 12, "I" * 16 + "5" * 12),
    ),
    "read2": fastq(
        ("p1", "GATTACAGATCGGAAGAGCGGTTCAGCAGG", "I" * 30),
        ("p2", "GATTACAGATCGGAAGAGCGGTTCAGCAGG", "I" * 30),
        ("p3", "C" * 30, "5" * 30),
    ),
}
result = generate_qc_report_tenx(tenx)
print(result)
```

- **Example output:**

```json
{
  "module": "sequencing_qc",
  "operation": "generate_qc_report_tenx",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 1,
    "text_length": 446,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "c2bd8879256656753e784d4d4fe4f7d2de6add20603c60c539095f494a8d364c"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "chemistry": "v3",
    "barcode_length": 16,
    "umi_length": 12,
    "read_pairs": 3,
    "sections": {
      "barcode": {
        "yield": {
          "...": "9 keys"
        },
        "basecalling": {
          "...": "12 keys"
        }
      },
      "umi": {
        "yield": {
          "...": "9 keys"
        },
        "basecalling": {
          "...": "12 keys"
        }
      },
      "rna_read": {
        "yield": {
          "...": "9 keys"
        },
        "basecalling": {
          "...": "12 keys"
        },
        "compact": {
          "...": "11 keys"
        }
      }
    },
    "barcode_and_umi": {
      "chemistry": "v3",
      "barcode_length": 16,
      "umi_length": 12,
      "read_pairs": 3,
      "cell_barcode": {
        "reads": 3,
        "barcode_source": "read_start",
        "distinct_barcodes": 2,
        "reads_with_N_in_barcode": 0,
        "top_barcodes": [
          {
            "...": "2 keys"
          },
          {
            "...": "2 keys"
          }
        ],
        "barcode_entropy_bits": 0.918296,
        "mean_quality": 40.0,
        "q30_fraction": 1.0,
        "mean_quality_by_position": [
          40.0,
          40.0,
          "..."
        ]
      },
      "umi": {
        "distinct_umis": 2,
        "reads_with_N_in_umi": 0,
        "mean_quality": 33.333333,
        "q30_fraction": 0.666667,
        "homopolymer_umi_reads": 3
      }
    },
    "duplication": {
      "distinct_sequences": 2,
      "definition": "Reads are duplicates when cell barcode, UMI and the first 25 bases of read 2 are identical. This is a proxy: Cell Ranger's sequencing saturation uses (barcode, ...",
      "duplicate_rate": 0.333333
    },
    "yield": {
      "chemistry": "v3",
      "barcode_length": 16,
      "umi_length": 12,
      "read_pairs": 3,
      "total_bases_read1_read2": 174,
      "distinct_barcodes": 2,
      "median_reads_per_barcode": 1.5
    },
    "report_text": "DNAKit 10x Chromium QC report\n=============================\nchemistry v3: 16 bp barcode + 12 bp UMI; read pairs: 3\nQ30 fraction (barcode / UMI / RNA read): 1.0 ...",
    "output_path": null
  }
}
```

## 44) plot_per_base_quality

- **Purpose:** Draws a per-base quality box plot (FastQC style) as an SVG string.

- **Method:** For every read position (positions are grouped when the reads are longer than `max_columns`) a box shows the lower to upper quartile of Phred Q, the red line the median, the whiskers the 10th to 90th percentile and the blue line the mean. The SVG is returned as text and written to a file only when `output_path` is given. Pure Python; no plotting library.

- **API:** `dnakit.sequencing_qc.plot_per_base_quality(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input` is FASTQ text or a file path (`.gz` allowed), a list of reads ((id, sequence, quality) tuples, `{sequence, quality}` mappings, FASTQ record objects or bare sequences), or a pair (then `mate` selects the read). Phred qualities are required. Options: `mate` (`read1` default or `read2`), `max_columns` (default 150), `output_path` (optional `.svg`; no overwrite unless `overwrite=True`), `phred_offset` (33 or 64, default 33), `max_reads` (default 1000000).

- **Output:** Returns a JSON-serializable envelope with `status` (`"ok"`), `module`, `operation`, `provenance` (`implementation: "native"`, `backend: "local"`), `input_summary`, `warnings` and `result`. `result` holds `svg` (the SVG text), `output_path` (`null` unless a file was written), `positions`, `positions_per_box` and a `legend`. Invalid input or options raise `InputFormatError` or `ConfigurationError`.

- **Example code:**

```python
from dnakit.sequencing_qc import plot_per_base_quality

def fastq(*records):
    return "".join(f"@{name}\n{seq}\n+\n{qual}\n" for name, seq, qual in records)

reads = fastq(
    ("r1", "ACGTACGTAC", "I" * 10),
    ("r2", "ACGTACGTAC", "I" * 10),
    ("r3", "TTTTTGGGGG", "5" * 10),
    ("r4", "NNACGTACGT", "!!" + "I" * 8),
)
result = plot_per_base_quality(reads)
print(result)
```

- **Example output:**

```json
{
  "module": "sequencing_qc",
  "operation": "plot_per_base_quality",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 4,
    "text_length": 112,
    "sequence_bases": 40,
    "lengths": [10, 10, 10, 10],
    "sha256": "74e1d5397e820db4ae7e5d14647765715a6c1e7a262dff425226e4bce9a7169a"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "svg": "<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"370\" height=\"378\" viewBox=\"0 0 370 378\" font-family=\"sans-serif\" font-size=\"11\">\n<title>Quality scores per base (...",
    "output_path": null,
    "positions": 10,
    "positions_per_box": 1,
    "legend": "box = lower to upper quartile, red line = median, whiskers = 10th-90th percentile, blue line = mean; background: green Q>=28, yellow 20-28, red <20."
  }
}
```

## 45) plot_qc_summary

- **Purpose:** Draws a bar chart of key QC percentages (Q20, Q30, GC, N, duplication, adapter) as an SVG string.

- **Method:** The bars are the percentage of bases at Q20 or above, at Q30 or above, GC content, N bases, duplicated reads (exact-sequence duplicates) and reads with an adapter hit (k-mer rule of `calculate_adapter_content`). The SVG is returned as text together with the metric values and written to a file only when `output_path` is given. Pure Python; no plotting library.

- **API:** `dnakit.sequencing_qc.plot_qc_summary(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input` is FASTQ text or a file path (`.gz` allowed), a list of reads ((id, sequence, quality) tuples, `{sequence, quality}` mappings, FASTQ record objects or bare sequences). Options: `adapters` (as `calculate_adapter_content`), `output_path` (optional `.svg`; no overwrite unless `overwrite=True`), `phred_offset` (33 or 64, default 33), `max_reads` (default 1000000).

- **Output:** Returns a JSON-serializable envelope with `status` (`"ok"`), `module`, `operation`, `provenance` (`implementation: "native"`, `backend: "local"`), `input_summary`, `warnings` and `result`. `result` holds `svg` (the SVG text), `output_path` (`null` unless a file was written) and `metrics` (percentages for `bases >= Q20`, `bases >= Q30`, `GC content`, `N bases`, `duplicated reads`, `reads with adapter`). Invalid input or options raise `InputFormatError` or `ConfigurationError`.

- **Example code:**

```python
from dnakit.sequencing_qc import plot_qc_summary

def fastq(*records):
    return "".join(f"@{name}\n{seq}\n+\n{qual}\n" for name, seq, qual in records)

reads = fastq(
    ("r1", "ACGTACGTAC", "I" * 10),
    ("r2", "ACGTACGTAC", "I" * 10),
    ("r3", "TTTTTGGGGG", "5" * 10),
    ("r4", "NNACGTACGT", "!!" + "I" * 8),
)
result = plot_qc_summary(reads)
print(result)
```

- **Example output:**

```json
{
  "module": "sequencing_qc",
  "operation": "plot_qc_summary",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 4,
    "text_length": 112,
    "sequence_bases": 40,
    "lengths": [10, 10, 10, 10],
    "sha256": "74e1d5397e820db4ae7e5d14647765715a6c1e7a262dff425226e4bce9a7169a"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "svg": "<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"660\" height=\"204\" viewBox=\"0 0 660 204\" font-family=\"sans-serif\" font-size=\"11\">\n<title>Sequencing QC summary</ti...",
    "output_path": null,
    "metrics": {
      "bases >= Q20": 95.0,
      "bases >= Q30": 70.0,
      "GC content": 50.0,
      "N bases": 5.0,
      "duplicated reads": 25.0,
      "reads with adapter": 0.0
    }
  }
}
```

## 46) plot_qc_summary_tenx

- **Purpose:** Draws a bar chart of 10x QC percentages (Q30 in barcode/UMI/RNA read, adapter, duplication proxy, GC) as an SVG string.

- **Method:** The bars are the Q30 fraction of the barcode, UMI and RNA read, the share of RNA reads with an adapter, template-switch-oligo or poly-A hit, the duplicate-rate proxy of `calculate_duplication_tenx`, the GC content of the RNA read and, with a `whitelist`, the fraction of valid barcodes. The SVG is returned as text together with the metric values and written to a file only when `output_path` is given. The 10x chemistry is not auto-detected: `chemistry` (or the explicit lengths) decides how read 1 is split into cell barcode and UMI.

- **API:** `dnakit.sequencing_qc.plot_qc_summary_tenx(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input` is a pair `{"read1": <barcode+UMI reads>, "read2": <cDNA reads>}` (FASTQ text, path or read list); a single read set is rejected and read 1 must be at least `barcode_length` + `umi_length` long. Options: `chemistry` (`v1` 14+10, `v2` 16+10, `v3` 16+12 bp barcode+UMI; default `v3`), `barcode_length`, `umi_length`, `phred_offset` (33), `max_reads` (1000000). Extra options: `whitelist` (adds the valid-barcode bar), `adapters`, `output_path` (optional `.svg`; no overwrite unless `overwrite=True`).

- **Output:** Returns a JSON-serializable envelope with `status` (`"ok"`), `module`, `operation`, `provenance` (`implementation: "native"`, `backend: "local"`), `input_summary`, `warnings` and `result`. `result` holds `chemistry`, `barcode_length`, `umi_length`, `read_pairs`, `svg` (the SVG text), `output_path` (`null` unless a file was written) and `metrics` (percentages such as `Q30 bases in barcode`, `Q30 bases in UMI`, `Q30 bases in RNA read`, `RNA reads with adapter/TSO/polyA`, `duplicate rate (proxy)`, `GC content of RNA read` and, with a whitelist, `valid barcodes`). Invalid input or options raise `InputFormatError` or `ConfigurationError`.

- **Example code:**

```python
from dnakit.sequencing_qc import plot_qc_summary_tenx

def fastq(*records):
    return "".join(f"@{name}\n{seq}\n+\n{qual}\n" for name, seq, qual in records)

barcode = "ACGTACGTACGTACGT"
tenx = {
    "read1": fastq(
        ("p1", barcode + "A" * 12, "I" * 28),
        ("p2", barcode + "A" * 12, "I" * 28),
        ("p3", "T" * 16 + "C" * 12, "I" * 16 + "5" * 12),
    ),
    "read2": fastq(
        ("p1", "GATTACAGATCGGAAGAGCGGTTCAGCAGG", "I" * 30),
        ("p2", "GATTACAGATCGGAAGAGCGGTTCAGCAGG", "I" * 30),
        ("p3", "C" * 30, "5" * 30),
    ),
}
result = plot_qc_summary_tenx(tenx, options={"whitelist": [barcode]})
print(result)
```

- **Example output:**

```json
{
  "module": "sequencing_qc",
  "operation": "plot_qc_summary_tenx",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 1,
    "text_length": 446,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "c2bd8879256656753e784d4d4fe4f7d2de6add20603c60c539095f494a8d364c"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "chemistry": "v3",
    "barcode_length": 16,
    "umi_length": 12,
    "read_pairs": 3,
    "svg": "<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"660\" height=\"228\" viewBox=\"0 0 660 228\" font-family=\"sans-serif\" font-size=\"11\">\n<title>10x Chromium QC summary (...",
    "output_path": null,
    "metrics": {
      "Q30 bases in barcode": 100.0,
      "Q30 bases in UMI": 66.6667,
      "Q30 bases in RNA read": 66.6667,
      "RNA reads with adapter/TSO/polyA": 66.6667,
      "duplicate rate (proxy)": 33.3333,
      "GC content of RNA read": 68.8889,
      "valid barcodes": 66.6667
    }
  }
}
```

## 47) run_cutadapt

- **Purpose:** Document a reusable run cutadapt capability for Sequencing quality control workflows. Candidate CAP-0683 comes from 1 source row(s) and is currently marked `new_candidate`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the local paths and typed options, invoke the caller-installed Cutadapt executable without a shell, enforce time/output limits, verify declared trimmed FASTQ files and an optional JSON report, and return structured provenance. The fixed option allowlist is a DNAKit boundary rather than a claim to expose the backend's complete CLI. The wrapper records the Martin (2011), DOI 10.14806/ej.17.1.200 citation and MIT license evidence.

- **API:** `dnakit.sequencing_qc.run_cutadapt(input, *, options=None, verbose=True)` (formal DNAKit namespace; the same bounded external-process effect is registered for Agent/MCP).

- **Input:** Required `input`: a FASTQ path or a `{read1, read2}` mapping. Required `options` values include the caller-installed `executable_path` and output location(s); `verbose` is optional. Only the documented allow-listed options are accepted, and output paths must not already exist.

- **Output:** The wrapper returns a JSON-serializable completion envelope with `status`, `backend`, the exact argv tokens, declared `outputs`, bounded backend output, and `provenance`; the backend produces declared trimmed FASTQ files and an optional JSON report. Missing executables, non-zero exits, missing artifacts, invalid JSON reports and output-limit violations are reported as structured errors.

- **Example code:**

```python
from pathlib import Path

from dnakit.sequencing_qc import run_cutadapt

result = run_cutadapt(
    input={"read1": Path("reads_R1.fastq"), "read2": Path("reads_R2.fastq")},
    options={
        "executable_path": "/usr/local/bin/cutadapt",
        "output1": Path("results/trimmed_R1.fastq"),
        "output2": Path("results/trimmed_R2.fastq"),
        "adapters": ["AGATCGGAAG"],
        "json_report": Path("results/cutadapt.json"),
    },
    verbose=False,
)
print(result)
```

- **Example output:**

```json
{
  "status": "completed",
  "backend": "cutadapt",
  "executable": "/usr/local/bin/backend",
  "arguments": [
    "<validated argv tokens>"
  ],
  "outputs": [
    "results/report.html"
  ],
  "provenance": {
    "implementation": "bounded_external_adapter",
    "citation": "Martin (2011), DOI 10.14806/ej.17.1.200",
    "license": "MIT"
  }
}
```

## 48) run_fastp

- **Purpose:** Document a reusable run fastp capability for Sequencing quality control workflows. Candidate CAP-0684 comes from 1 source row(s) and is currently marked `new_candidate`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the local paths and typed options, invoke the caller-installed fastp executable without a shell, enforce time/output limits, verify declared filtered FASTQ files plus JSON and HTML reports, and return structured provenance. The fixed option allowlist is a DNAKit boundary rather than a claim to expose the backend's complete CLI. The wrapper records the Chen (2025), DOI 10.1002/imt2.70078 citation and MIT license evidence.

- **API:** `dnakit.sequencing_qc.run_fastp(input, *, options=None, verbose=True)` (formal DNAKit namespace; the same bounded external-process effect is registered for Agent/MCP).

- **Input:** Required `input`: a FASTQ path or a `{read1, read2}` mapping. Required `options` values include the caller-installed `executable_path` and output location(s); `verbose` is optional. Only the documented allow-listed options are accepted, and output paths must not already exist.

- **Output:** The wrapper returns a JSON-serializable completion envelope with `status`, `backend`, the exact argv tokens, declared `outputs`, bounded backend output, and `provenance`; the backend produces declared filtered FASTQ files plus JSON and HTML reports. Missing executables, non-zero exits, missing artifacts, invalid JSON reports and output-limit violations are reported as structured errors.

- **Example code:**

```python
from pathlib import Path

from dnakit.sequencing_qc import run_fastp

result = run_fastp(
    input=Path("reads.fastq"),
    options={
        "executable_path": "/usr/local/bin/fastp",
        "output": Path("results/filtered.fastq"),
        "json_report": Path("results/fastp.json"),
        "html_report": Path("results/fastp.html"),
    },
    verbose=False,
)
print(result)
```

- **Example output:**

```json
{
  "status": "completed",
  "backend": "fastp",
  "executable": "/usr/local/bin/backend",
  "arguments": [
    "<validated argv tokens>"
  ],
  "outputs": [
    "results/report.html"
  ],
  "provenance": {
    "implementation": "bounded_external_adapter",
    "citation": "Chen (2025), DOI 10.1002/imt2.70078",
    "license": "MIT"
  }
}
```

## 49) run_fastqc

- **Purpose:** Document a reusable run fastqc capability for Sequencing quality control workflows. Candidate CAP-0685 comes from 1 source row(s) and is currently marked `new_candidate`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the local paths and typed options, invoke the caller-installed FastQC 0.12.x executable without a shell, enforce time/output limits, verify FastQC HTML and ZIP reports in an existing output directory, and return structured provenance. The fixed option allowlist is a DNAKit boundary rather than a claim to expose the backend's complete CLI. The wrapper records the FastQC project documentation; no DOI supplied by the survey citation and GPL-3.0 license evidence.

- **API:** `dnakit.sequencing_qc.run_fastqc(input, *, options=None, verbose=True)` (formal DNAKit namespace; the same bounded external-process effect is registered for Agent/MCP).

- **Input:** Required `input`: one to eight local FASTQ/BAM/SAM-style read paths. Required `options` values include the caller-installed `executable_path` and output location(s); `verbose` is optional. Only the documented allow-listed options are accepted, and output paths must not already exist.

- **Output:** The wrapper returns a JSON-serializable completion envelope with `status`, `backend`, the exact argv tokens, declared `outputs`, bounded backend output, and `provenance`; the backend produces FastQC HTML and ZIP reports in an existing output directory. Missing executables, non-zero exits, missing artifacts, invalid JSON reports and output-limit violations are reported as structured errors.

- **Example code:**

```python
from pathlib import Path

from dnakit.sequencing_qc import run_fastqc

result = run_fastqc(
    input=Path("reads.fastq.gz"),
    options={
        "executable_path": "/opt/FastQC/fastqc",
        "output_directory": Path("results/fastqc"),
        "kmer": 7,
    },
    verbose=False,
)
print(result)
```

- **Example output:**

```json
{
  "status": "completed",
  "backend": "fastqc",
  "executable": "/usr/local/bin/backend",
  "arguments": [
    "<validated argv tokens>"
  ],
  "outputs": [
    "results/report.html"
  ],
  "provenance": {
    "implementation": "bounded_external_adapter",
    "citation": "FastQC project documentation; no DOI supplied by the survey",
    "license": "GPL-3.0"
  }
}
```

## 50) run_multiqc

- **Purpose:** Document a reusable run multiqc capability for Sequencing quality control workflows. Candidate CAP-0686 comes from 1 source row(s) and is currently marked `new_candidate`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the local paths and typed options, invoke the caller-installed MultiQC executable without a shell, enforce time/output limits, verify a declared MultiQC HTML report in an empty output directory, and return structured provenance. The fixed option allowlist is a DNAKit boundary rather than a claim to expose the backend's complete CLI. The wrapper records the Ewels et al. (2016), DOI 10.1093/bioinformatics/btw354 citation and GPL-3.0-or-later license evidence.

- **API:** `dnakit.sequencing_qc.run_multiqc(input, *, options=None, verbose=True)` (formal DNAKit namespace; the same bounded external-process effect is registered for Agent/MCP).

- **Input:** Required `input`: one or more local files/directories containing QC inputs. Required `options` values include the caller-installed `executable_path` and output location(s); `verbose` is optional. Only the documented allow-listed options are accepted, and output paths must not already exist.

- **Output:** The wrapper returns a JSON-serializable completion envelope with `status`, `backend`, the exact argv tokens, declared `outputs`, bounded backend output, and `provenance`; the backend produces a declared MultiQC HTML report in an empty output directory. Missing executables, non-zero exits, missing artifacts, invalid JSON reports and output-limit violations are reported as structured errors.

- **Example code:**

```python
from pathlib import Path

from dnakit.sequencing_qc import run_multiqc

result = run_multiqc(
    input=[Path("results/fastqc")],
    options={
        "executable_path": "/opt/multiqc/bin/multiqc",
        "output_directory": Path("results/multiqc"),
        "filename": "multiqc_report.html",
    },
    verbose=False,
)
print(result)
```

- **Example output:**

```json
{
  "status": "completed",
  "backend": "multiqc",
  "executable": "/usr/local/bin/backend",
  "arguments": [
    "<validated argv tokens>"
  ],
  "outputs": [
    "results/report.html"
  ],
  "provenance": {
    "implementation": "bounded_external_adapter",
    "citation": "Ewels et al. (2016), DOI 10.1093/bioinformatics/btw354",
    "license": "GPL-3.0-or-later"
  }
}
```
