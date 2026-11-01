# khmer bounded read parsing

## khmer_parse_reads

**Purpose：** Read FASTA/FASTQ records, preserving names, sequences and FASTQ qualities.

**Method：** Calls the actual khmer.ReadParser; pairing uses iter_read_pairs. Limits stop materialization and set truncated.

**API：** `dnakit.microbiome.khmer_p27_read_parser_functions.khmer_parse_reads(path: 'str | Path', *, max_records: 'int' = 1000, max_bases: 'int' = 1000000, verbose: bool = True) -> 'dict[str, Any]'`

**Input：** `path`, `max_records=1000`, `max_bases=1000000`, `verbose=True`.

**Output：** `records`, `record_count`, `base_count`, `truncated`, `provenance`.

**Example file:** `pairs.fa`

```text
>read/1
ACGT
>read/2
TGCA
```

**Example code:**

```python
from dnakit.microbiome.khmer_p27_read_parser_functions import khmer_parse_reads
result = khmer_parse_reads("pairs.fa", verbose=False)
print(result["record_count"], result["base_count"], result["truncated"])
```

**Example output:**

```text
2 8 False
```

## khmer_parse_read_pairs

**Purpose：** Read interleaved FASTA/FASTQ pairs with khmer's name-pairing validation.

**Method：** Calls the actual khmer.ReadParser; pairing uses iter_read_pairs. Limits stop materialization and set truncated.

**API：** `dnakit.microbiome.khmer_p27_read_parser_functions.khmer_parse_read_pairs(path: 'str | Path', *, max_records: 'int' = 1000, max_bases: 'int' = 1000000, verbose: bool = True) -> 'dict[str, Any]'`

**Input：** `path`, `max_records=1000`, `max_bases=1000000`, `verbose=True`.

**Output：** `records`, `record_count`, `base_count`, `truncated`, `provenance`.

**Example file:** `pairs.fa`

```text
>read/1
ACGT
>read/2
TGCA
```

**Example code:**

```python
from dnakit.microbiome.khmer_p27_read_parser_functions import khmer_parse_read_pairs
result = khmer_parse_read_pairs("pairs.fa", verbose=False)
print(result["record_count"], result["base_count"], result["truncated"])
```

**Example output:**

```text
1 8 False
```

