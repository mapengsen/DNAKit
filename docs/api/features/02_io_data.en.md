# File reading and writing

Read and write common DNA sequences, annotations, tables, and compressed files. Ordinary users only need to remember `dnakit.read()` when reading a sequence file, and only need to remember `dnakit.write()` when writing out; when the file is large, the same method is still used, and only the reading parameter is set to `mode="stream"`.

Large files are not a new file format and do not require new common read and write entries. The original `IO-005` number continues to be used for requirements tracking, and its streaming, chunking, and indexing capabilities are unified into the large file usage of `IO-001`.

The read interface preserves input order by default.

## 1) IO-001 sequence format

- **Function:** Read and write out FASTA, FASTQ, GenBank files, convert sequences and available IDs, quality values and annotations into DNAKit records for continued analysis or format conversion.
- **Unified API:** `dnakit.read(source[required], format[optional], config[optional], mode[optional])`, `dnakit.write(records[required], target[required], format[optional], config[optional])`; set `mode="dna"` for normal reading, `mode="stream"` for large files.
- **Compatible API:** `read_one()`, `read_set()` are still available, but new code does not need to remember these names.
- **Input:** `mode="dna"` Returns `DNA` whether the file contains one or multiple records; `mode="stream"` returns a single consumption data stream. Writing can be passed directly to `DNA`, or streaming records or other record iterators can be passed directly. There is no need for a special "large file writing" method.
- **Sample code:**

```python
from io import StringIO

from dnakit import read, write

source = StringIO(">seq1 demo\nACGTN\n")
records = read(source, format="fasta", mode="dna")

output = StringIO()
written = write(records, output, format="fasta")
print(type(records).__name__, records.ids, records[0].symbols)
print(written.record_count)
```

- **Example results:**

```text
DNA ('seq1',) ACGTN
1
```

### Large file parameters and advanced indexing (Track: IO-005) {#5-io-005}

- **Ordinary users:** Use the same `dnakit.read(..., mode="stream")` to read one by one, and then pass the returned data stream directly to `dnakit.write()`; use `iter_chunks()` only when fixed chunking is required.

```python
import dnakit

with dnakit.read("large.fa", mode="stream") as records:
    dnakit.write(records, "copy.fa")
```

- **Advanced users:** Only use the index interface when you need to randomly access large FASTA/FASTQ by record ID or coordinates; these interfaces are not required for ordinary file reading.
- **Advanced API:** `dnakit.io.iter_chunks(values[required], chunk_size[required])`, `dnakit.io.build_fasta_index(source_path[required], index_path[optional], overwrite[optional], max_records[optional], max_line_length[optional])`, `dnakit.io.load_fasta_index(index_path[required], source_path[optional], verify_checksum[optional], max_index_bytes[optional])`, `dnakit.io.build_fastq_index(source_path[required], index_path[optional], overwrite[optional], phred_offset[optional], max_records[optional], max_line_length[optional], max_record_bytes[optional], max_source_bytes[optional])`, `dnakit.io.load_fastq_index(index_path[required], source_path[optional], verify_checksum[optional], max_index_bytes[optional], max_entries[optional])`, `dnakit.io.FastaIndex.fetch(record_id[required], start[optional], end[optional], strand[optional], max_record_bytes[optional])`, `dnakit.io.FastqIndex.fetch(record_id[required], start[optional], end[optional], strand[optional], max_record_bytes[optional])`.

## 2) IO-002 comment format

- **Function:** Read and write GFF3, BED, and AGP annotation files, and unify the coordinate rules of different formats to generate structured data that can be used for region query, sequence segmentation, and annotation conversion.
- **API:** `dnakit.io.read_gff3(source[required], max_records[optional], max_line_length[optional], max_header_lines[optional])`, `dnakit.io.write_gff3(document[required], target[required], overwrite[optional], create_parents[optional])`, `dnakit.io.read_bed(source[required], max_records[optional], max_line_length[optional], max_header_lines[optional])`, `dnakit.io.write_bed(document[required], target[required], overwrite[optional], create_parents[optional])`, `dnakit.io.read_agp(source[required], max_records[optional], max_line_length[optional], max_header_lines[optional])`, `dnakit.io.write_agp(document[required], target[required], overwrite[optional], create_parents[optional])`.
- **Input:** The annotation file path or text stream is required when reading, and the number of records, line length and header limit are optional; when writing, pass in the corresponding document/entry and target path or text stream.
- **Sample code:**

```python
from io import StringIO

from dnakit.io import read_gff3, write_gff3

source = StringIO(
    "##gff-version 3\n"
    "chr1\ttest\tgene\t2\t4\t.\t+\t.\tID=g1\n"
)
document = read_gff3(source)

output = StringIO()
count = write_gff3(document, output)
print(document.entries[0].feature.location)
print(count)
```

- **Example results:**

```text
Interval(start=1, end=4)
1
```

## 3) IO-003 table format

- **Function:** Read and write CSV, TSV, JSON and Parquet according to the specified table structure, and perform reproducible conversion between DNA records, metadata and general tables.
- **API:** `dnakit.io.TableSchema(columns[required], schema_version[optional], column_types[optional], nullable[optional])`, `dnakit.io.export_table(rows[required], target[required], format[required], schema[required], overwrite[optional], max_rows[optional], max_columns[optional], max_cell_characters[optional], max_output_bytes[optional], null_value[optional], parquet_compression[optional])`, `dnakit.io.read_table(source[required], format[required], schema[required], missing_values[optional], max_rows[optional], max_columns[optional], max_cell_characters[optional], max_file_bytes[optional], max_decoded_bytes[optional])`, `dnakit.io.export_result(result[required], target[required], format[optional], overwrite[optional])`, `dnakit.io.parquet_backend_status()`.
- **Input:** `export_table()`/`read_table()` uses 2D row data or target path, explicit `format` and `TableSchema`, optional column types, nullable fields, missing value flags, compression method and resource cap; `export_result()` receives a single result object with `to_dict()` and target path, defaults to JSON and automatically derives a single row schema.
- **Sample code:**

```python
from pathlib import Path
from tempfile import TemporaryDirectory

from dnakit.io import TableSchema, export_table, read_table

schema = TableSchema(
    ("id", "score"),
    column_types={"id": "string", "score": "number"},
    nullable=(),
)
with TemporaryDirectory() as directory:
    path = Path(directory) / "scores.json"
    exported = export_table(
        [{"id": "a", "score": 0.75}],
        path,
        format="json",
        schema=schema,
    )
    loaded = read_table(path, format="json", schema=schema)
    print(exported.row_count, loaded.rows[0]["score"])
```

- **Example results:**

```text
1 0.75
```

## 4) IO-004 compressed file

- **Function:** Directly read and write gzip-compressed sequence or table files without the caller decompressing them first, and retain the same parsing results as ordinary files.
- **API:** `dnakit.read(source[required], format[optional], config[optional], mode[optional])`, `dnakit.write(records[required], target[required], format[optional], config[optional])`; the old `read_one()` is only a compatible entry.
- **Input:** Required compressed file path or binary-enabled stream; use `compression="auto"|"none"|"gzip"` to control compression, optionally writing out the compression level.
- **Sample code:**

```python
from pathlib import Path
from tempfile import TemporaryDirectory

from dnakit import DNA, ReadConfig, WriteConfig, read, write

with TemporaryDirectory() as directory:
    path = Path(directory) / "seq.fa.gz"
    record = DNA("ACGT", id="seq1")
    write(
        record,
        path,
        format="fasta",
        config=WriteConfig(compression="gzip"),
    )
    loaded = read(
        path,
        format="fasta",
        config=ReadConfig(compression="auto"),
        mode="dna",
    )
    print(loaded.id, loaded.symbols)
```

- **Example results:**

```text
seq1 ACGT
```
