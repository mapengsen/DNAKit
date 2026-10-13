# G1b: offline pyfaidx, dnaio and screed wrappers

This group selects the DNA, FASTA/FASTQ and interval-parsing functions from the public inventories of pyfastx, pyfaidx, dnaio, screed and py2bit, and wraps them as offline functions in `dnakit.sequence`. Eight new wrappers are added: three from pyfaidx, two from dnaio and three from screed. pyfastx and py2bit need no new wrapper. Every other inventory entry has a recorded reason, such as a duplicate of an existing dnakit function, an existing wrapper in `sequence/backends.py`, a command-line entry point, an internal helper, a deprecated alias or upstream test code. All functions read only explicit local files or strings and never use the network.

## Coverage

The table is counted from the ledger; statuses are defined there.

| Package | New wrappers | Already wrapped | Duplicate of dnakit | Family variant | Skipped |
| --- | --- | --- | --- | --- | --- |
| pyfastx | 0 | 0 | 1 | 0 | 2 |
| pyfaidx | 3 | 0 | 2 | 0 | 15 |
| dnaio | 2 | 0 | 1 | 0 | 6 |
| screed | 3 | 3 | 16 | 5 | 68 |
| py2bit | 0 | 1 | 0 | 0 | 0 |

## Example data

The examples import the functions from the package. The files only demonstrate the interface; paired read names follow the `/1` and `/2` convention. The screed database is built in a temporary directory.

**Example file:** `reads.fa`

```text
>r1 first read
ACGTNACGT
>r2 second read
ggcc
```

**Example file:** `reads.fq`

```text
@q1 desc
ACGTAC
+
IIII#!
@q2
GGCCAA
+
ABCDEF
```

**Example file:** `mate_1.fq`

```text
@pair7/1 lane1
ACGTAC
+
IIIIII
@pair8/1
GGCC
+
IIII
```

**Example file:** `mate_2.fq`

```text
@pair7/2 lane1
TTGCAA
+
IIIIII
@pair8/2
CCAA
+
IIII
```

## 1) Translate chromosome names with explicit pairs

- **Purpose:** Translate chromosome names through an explicit old-to-new table.
- **Method:** Calls pyfaidx.translate_chr_name to build the mapping. The two name lists must have the same length, and old names must not repeat. Names without a pair go to `unmapped`; nothing is guessed.
- **API:** `dnakit.sequence.translate_chromosome_names(names[required], *, from_names[required], to_names[required], verbose[optional])`
- **Input:** A list of names to translate, plus two equal-length lists of old and new names.
- **Output:** `translations` lists `name` and `translated` in input order; `unmapped` is the sorted list of unmatched names; `requested_count` and `translated_count` are included.
- **Example code:**

```python
from dnakit.sequence import translate_chromosome_names

result = translate_chromosome_names(
    ["chr1", "chrM", "chr2", "chrX"],
    from_names=["chr1", "chr2", "chrM"],
    to_names=["1", "2", "MT"],
    verbose=False,
)
print(result["translations"])
print(result["unmapped"])
```

**Example output:**

```text
[{'name': 'chr1', 'translated': '1'}, {'name': 'chrM', 'translated': 'MT'}, {'name': 'chr2', 'translated': '2'}]
['chrX']
```

## 2) Split one BED interval line

- **Purpose:** Parse one BED line into chromosome, 0-based start and half-open end.
- **Method:** Calls pyfaidx.bed_split. Lines starting with `#` are comments. Negative starts, ends before starts, blank lines and malformed lines raise an error.
- **API:** `dnakit.sequence.split_bed_entry(line[required], verbose[optional])`
- **Input:** One BED line with at least three columns: chromosome, start and end.
- **Output:** `is_comment`, `chromosome`, `start`, `end` and `length`. Coordinate fields are `None` for comment lines.
- **Example code:**

```python
from dnakit.sequence import split_bed_entry

entry = split_bed_entry("chr1\t10\t20\tgeneA", verbose=False)
print(entry["chromosome"], entry["start"], entry["end"], entry["length"])
print(split_bed_entry("# header", verbose=False)["is_comment"])
```

**Example output:**

```text
chr1 10 20 10
True
```

## 3) Parse a UCSC region

- **Purpose:** Convert a UCSC `chr:start-end` region into a 0-based half-open interval; a bare name means the whole sequence.
- **Method:** Calls pyfaidx.ucsc_split and subtracts one from the 1-based inclusive start. Unparseable coordinates, such as thousands separators, raise an error instead of being read as a whole sequence. Names that contain a colon need an explicit interval.
- **API:** `dnakit.sequence.split_ucsc_region(region[required], verbose[optional])`
- **Input:** One region string, for example `chr1:101-200` or `chr7`.
- **Output:** `chromosome`, `start`, `end`, `length` and `is_whole_chromosome`. Coordinates are `None` for a whole sequence.
- **Example code:**

```python
from dnakit.sequence import split_ucsc_region

region = split_ucsc_region("chr1:101-200", verbose=False)
print(region["chromosome"], region["start"], region["end"], region["length"])
print(split_ucsc_region("chr7", verbose=False)["is_whole_chromosome"])
```

**Example output:**

```text
chr1 100 200 100
True
```

## 4) Read FASTA/FASTQ and paired reads

- **Purpose:** Read FASTA or FASTQ with dnaio, including two-file paired FASTQ input.
- **Method:** Calls `dnaio.open`. `max_records` bounds the count, and pairs are counted as one entry. dnaio rejects mates whose names disagree; the wrapper reports that as ConfigurationError. Use `interleaved=True` for one interleaved file.
- **API:** `dnakit.sequence.read_fastx_records(source[required], *, mate_source[optional], interleaved[optional], max_records[optional], verbose[optional])`
- **Input:** A FASTA or FASTQ path; optionally a second mate file in `mate_source`, or `interleaved=True` for one interleaved file.
- **Output:** `layout`, `format`, `entry_unit`, `entry_count` and `records`. Single-end records carry `name`, `id`, `comment`, `sequence` and `length`, plus `qualities` for FASTQ; paired records carry `read1` and `read2`.
- **Example code:**

```python
from dnakit.sequence import read_fastx_records

result = read_fastx_records("reads.fq", verbose=False)
print(result["format"], result["entry_count"])
for record in result["records"]:
    print(record["id"], record["sequence"], record["qualities"])
paired = read_fastx_records("mate_1.fq", mate_source="mate_2.fq", verbose=False)
print(paired["layout"], paired["entry_count"])
for pair in paired["records"]:
    print(pair["read1"]["id"], pair["read2"]["id"], pair["read1"]["length"])
```

**Example output:**

```text
fastq 2
q1 ACGTAC IIII#!
q2 GGCCAA ABCDEF
paired-end 2
pair7/1 pair7/2 6
pair8/1 pair8/2 4
```

## 5) Check paired read names

- **Purpose:** Decide whether a set of read names are mates of each other.
- **Method:** Compares names only: each name becomes a dnaio.SequenceRecord with no sequence or qualities, then `dnaio.records_are_mates` is called. Empty names and names with line breaks raise an error.
- **API:** `dnakit.sequence.check_read_mates(names[required], *, verbose[optional])`
- **Input:** A list of at least two read names.
- **Output:** `mates` as a boolean, plus `read_count` and the `names` as given.
- **Example code:**

```python
from dnakit.sequence import check_read_mates

print(check_read_mates(["pair7/1", "pair7/2"], verbose=False)["mates"])
print(check_read_mates(["pair7/1", "pair8/2"], verbose=False)["mates"])
```

**Example output:**

```text
True
False
```

## 6) Export a screed database

- **Purpose:** Write every record of a screed SQLite database to a FASTA or FASTQ file.
- **Method:** Counts records first: an empty database or one above `max_records` raises an error. Then calls `screed.ToFasta` or `screed.ToFastq`. An existing output file is never overwritten, and partial output is removed on failure or a count mismatch. screed wraps sequence lines at 80 columns.
- **API:** `dnakit.sequence.export_screed_database(database[required], output[required], *, output_format[optional], max_records[optional], verbose[optional])`
- **Input:** The database path returned by `build_screed_database` and a new output path. `output_format` is `fasta` (default) or `fastq`.
- **Output:** `record_count`, `output_format`, `database` and `output`.
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.sequence import build_screed_database
from dnakit.sequence import export_screed_database

work = Path(tempfile.mkdtemp())
built = build_screed_database("reads.fa", database=work / "reads", verbose=False)
exported = export_screed_database(built["database"], work / "out.fa", verbose=False)
print(exported["record_count"], exported["output_format"])
print((work / "out.fa").read_text(encoding="utf-8"), end="")
```

**Example output:**

```text
2 fasta
>r1 first read
ACGTNACGT
>r2 second read
ggcc
```

## 7) Parse FASTA/FASTQ with screed iterators

- **Purpose:** Parse an uncompressed FASTA or FASTQ file with screed's own record iterators.
- **Method:** `file_format` selects `screed.fasta.fasta_iter` or `screed.fastq.fastq_iter`. screed stores FASTQ header text as `annotations`; the wrapper reports it as `description`. A format mismatch or more than `max_records` records raises an error.
- **API:** `dnakit.sequence.read_screed_fastx_records(source[required], *, file_format[optional], parse_description[optional], max_records[optional], verbose[optional])`
- **Input:** An uncompressed FASTA or FASTQ file, with its format declared in `file_format`.
- **Output:** `format`, `record_count` and `records`. Each record has `name`, `description`, `length` and `sequence`; FASTQ records also have `quality`.
- **Example code:**

```python
from dnakit.sequence import read_screed_fastx_records

result = read_screed_fastx_records("reads.fa", verbose=False)
print(result["format"], result["record_count"])
for record in result["records"]:
    print(record["name"], record["description"], record["length"])
fastq = read_screed_fastx_records("reads.fq", file_format="fastq", verbose=False)
print(fastq["records"][0]["quality"])
```

**Example output:**

```text
fasta 2
r1 first read 9
r2 second read 4
IIII#!
```

## 8) Format records as FASTA/FASTQ text

- **Purpose:** Format explicit records as FASTA or FASTQ text without writing a file.
- **Method:** Each record is written to an in-memory buffer by `screed.screedRecord.write_fastx`. Names cannot contain whitespace; sequences may contain only ASCII letters, `-` and `*`; FASTQ qualities must match the sequence length and be printable ASCII. Sequences are not wrapped.
- **API:** `dnakit.sequence.format_fastx_records(records[required], *, output_format[optional], max_records[optional], verbose[optional])`
- **Input:** A list of dictionaries with `name` and `sequence`, optional `description`, and `quality` for `fastq` output.
- **Output:** `text`, `record_count` and `output_format`.
- **Example code:**

```python
from dnakit.sequence import format_fastx_records

fasta = format_fastx_records(
    [
        {"name": "r1", "description": "first read", "sequence": "ACGTN"},
        {"name": "r2", "sequence": "ggcc"},
    ],
    verbose=False,
)
print(fasta["text"], end="")
fastq = format_fastx_records(
    [{"name": "q1", "description": "desc", "sequence": "ACGTAC", "quality": "IIII#!"}],
    output_format="fastq",
    verbose=False,
)
print(fastq["text"], end="")
```

**Example output:**

```text
>r1 first read
ACGTN
>r2
ggcc
@q1 desc
ACGTAC
+
IIII#!
```

## Notes

- The wrappers are exported from `dnakit.sequence`.
- Upstream `screed.screedRecord.write_fastx_pair` calls an undefined name, `write_record`, and raises `NameError`, so it is not wrapped.
- On an empty screed database, upstream `ToFasta` and `ToFastq` raise `TypeError`. `export_screed_database` counts records first and reports a clear error.
- `pyfaidx.ucsc_split` returns a negative start for `chr1:0-5`; `split_ucsc_region` reports that as an error.
- `pyfastx.reverse_complement` does not validate its input, whereas `reverse_complement_dna` rejects non-IUPAC symbols.
