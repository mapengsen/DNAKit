# File reading and writing

Read with `dnakit.read()` and write with `dnakit.write()`.

## 1) sequence format {#1-io-001-sequence-format}

- **Purpose:** Read and write out FASTA, FASTQ, GenBank files, convert sequences and available IDs, quality values and annotations into DNAKit records for continued analysis or format conversion.
- **Method:** Parse records in file order; mode="stream" reads lazily and writing serializes the target format.
- **API:** `dnakit.read(source[required], format[optional], config[optional], mode[optional])`, `dnakit.write(records[required], target[required], format[optional], config[optional])`; set `mode="dna"` for normal reading, `mode="stream"` for large files.
- **Compatible API:** `read_one()`, `read_set()` are still available, but new code does not need to remember these names.
- **Input:** `mode="dna"` Returns `DNA` whether the file contains one or multiple records; `mode="stream"` returns a single consumption data stream. Writing can be passed directly to `DNA`, or streaming records or other record iterators can be passed directly. There is no need for a special "large file writing" method.
- **Output:** `read` returns `DNA` or a streaming `RecordSource`; `write` returns a `WriteResult` with the record count.
- **Example code:**

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

- **Example output:**

```text
DNA ('seq1',) ACGTN
1
```

<span id="5-io-005"></span>[FAQ: details](../../faq/tutorial_basics.md#large-files)

## 2) annotation formats {#2-io-002-annotation-formats}

- **Purpose:** Read and write GFF3, BED, and AGP annotation files, and unify the coordinate rules of different formats to generate structured data that can be used for region query, sequence segmentation, and annotation conversion.
- **Method:** Parse GFF3/BED/AGP fields and convert their coordinates into internal intervals.
- **API:** `dnakit.io.read_gff3(source[required], max_records[optional], max_line_length[optional], max_header_lines[optional])`, `dnakit.io.write_gff3(document[required], target[required], overwrite[optional], create_parents[optional])`, `dnakit.io.read_bed(source[required], max_records[optional], max_line_length[optional], max_header_lines[optional])`, `dnakit.io.write_bed(document[required], target[required], overwrite[optional], create_parents[optional])`, `dnakit.io.read_agp(source[required], max_records[optional], max_line_length[optional], max_header_lines[optional])`, `dnakit.io.write_agp(document[required], target[required], overwrite[optional], create_parents[optional])`.
- **Input:** The annotation file path or text stream is required when reading, and the number of records, line length and header limit are optional; when writing, pass in the corresponding document/entry and target path or text stream.
- **Output:** `read_gff3` returns an `AnnotationDocument`; `write_gff3` returns the written feature count.
- **Example code:**

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

- **Example output:**

```text
Interval(start=1, end=4)
1
```

## 3) table format {#3-io-003-table-format}

- **Purpose:** Read and write CSV, TSV, JSON and Parquet according to the specified table structure, and perform reproducible conversion between DNA records, metadata and general tables.
- **Method:** Validate columns, types and missing values against an explicit schema before parsing or serialization.
- **API:** `dnakit.io.TableSchema(columns[required], schema_version[optional], column_types[optional], nullable[optional])`, `dnakit.io.export_table(rows[required], target[required], format[required], schema[required], overwrite[optional], max_rows[optional], max_columns[optional], max_cell_characters[optional], max_output_bytes[optional], null_value[optional], parquet_compression[optional])`, `dnakit.io.read_table(source[required], format[required], schema[required], missing_values[optional], max_rows[optional], max_columns[optional], max_cell_characters[optional], max_file_bytes[optional], max_decoded_bytes[optional])`, `dnakit.io.export_result(result[required], target[required], format[optional], overwrite[optional])`, `dnakit.io.parquet_backend_status()`.
- **Input:** `export_table()`/`read_table()` uses 2D row data or target path, explicit `format` and `TableSchema`, optional column types, nullable fields, missing value flags, compression method and resource cap; `export_result()` receives a single result object with `to_dict()` and target path, defaults to JSON and automatically derives a single row schema.
- **Output:** `export_table` returns the written row count; `read_table` returns table rows and schema information.
- **Example code:**

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

- **Example output:**

```text
1 0.75
```

## 4) compressed file {#4-io-004-compressed-file}

- **Purpose:** Directly read and write gzip-compressed sequence or table files without the caller decompressing them first, and retain the same parsing results as ordinary files.
- **Method:** Decode or encode gzip around the data stream while reusing the normal format reader/writer.
- **API:** `dnakit.read(source[required], format[optional], config[optional], mode[optional])`, `dnakit.write(records[required], target[required], format[optional], config[optional])`; the old `read_one()` is only a compatible entry.
- **Input:** Required compressed file path or binary-enabled stream; use `compression="auto"|"none"|"gzip"` to control compression, optionally writing out the compression level.
- **Output:** Reading returns `DNA` or `RecordSource`; writing returns a file result using the selected compression.
- **Example code:**

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

- **Example output:**

```text
seq1 ACGT
```

## 5) Read and write alignment files {#bio-2}

- **Purpose:** Convert alignments between FASTA, Clustal, Stockholm and other formats.
- **Method:** Decode and encode one alignment with Biopython while preserving IDs.
- **API:** `dnakit.bioinformatics.alignment.read_alignment(path, format='fasta')`; `dnakit.bioinformatics.alignment.write_alignment(alignment, path, format='fasta')`
- **Input:** An alignment object or existing file; requires `bio-formats`. This example creates its input.
- **Output:** Reading returns `MultipleAlignment`; writing creates an alignment file in the selected format.
- **Example code:**

```python
from tempfile import TemporaryDirectory
from pathlib import Path
from dnakit.bioinformatics.alignment import MultipleAlignment, read_alignment, write_alignment

with TemporaryDirectory() as folder:
    msa = MultipleAlignment.from_sequences({"a": "ACGT", "b": "AC-T"})
    path = write_alignment(msa, Path(folder) / "aligned.fa")
    restored = read_alignment(path)
    print(restored.ids, restored.sequences)
```

- **Example output:**

```text
('a', 'b') ('ACGT', 'AC-T')
```

[Dependencies and implementation notes](../../faq/bioinformatics.md#alignment).

## 6) Read and write tree files {#bio-7-io}

- **Purpose:** Save a tree in Newick or another supported format and read its structure.
- **Method:** Use Biopython to decode and encode branch lengths and leaf nodes.
- **API:** `dnakit.bioinformatics.phylogeny.read_tree(path, format='newick')`; `dnakit.bioinformatics.phylogeny.write_tree(tree, path, format='newick')`
- **Input:** A tree object or existing tree file; requires `bio-formats`. This example creates its own tree.
- **Output:** Reading returns a `TreeNode` tree; writing creates a tree file in the selected format.
- **Example code:**

```python
from tempfile import TemporaryDirectory
from pathlib import Path
from dnakit.bioinformatics.phylogeny import TreeNode, read_tree, write_tree

tree = TreeNode(children=(TreeNode("a", 1), TreeNode("b", 1)))
with TemporaryDirectory() as folder:
    path = write_tree(tree, Path(folder) / "tree.nwk")
    restored = read_tree(path)
    print(restored.tips)
```

- **Example output:**

```text
('a', 'b')
```

[Dependencies and implementation notes](../../faq/bioinformatics.md#alignment).

## 7) Read FASTA/FASTQ and convert quality encoding {#bio-9}

- **Purpose:** Read sequencing records, save sequences and explicitly convert Phred+64/33.
- **Method:** Stream FASTA/FASTQ, including gzip files, with integer Phred qualities.
- **API:** `dnakit.bioinformatics.reads.Read(id, sequence, qualities=None)`; `dnakit.bioinformatics.reads.read_fastx(path, format='fastq', phred_offset=33)`; `dnakit.bioinformatics.reads.write_fastx(reads, path, format='fastq', phred_offset=33)`
- **Input:** File paths or `Read` records; FASTA-to-FASTQ conversion requires explicit `default_quality`.
- **Output:** Reading returns an iterator of `Read` records; conversion returns records with updated quality scores.
- **Example code:**

```python
from tempfile import TemporaryDirectory
from pathlib import Path
from dnakit.bioinformatics.reads import Read, read_fastx, write_fastx

with TemporaryDirectory() as folder:
    source = Path(folder) / "phred64.fq"
    target = Path(folder) / "phred33.fq"
    write_fastx([Read("r1", "ACGT", (30, 31, 32, 33))], source, phred_offset=64)
    count = write_fastx(read_fastx(source, phred_offset=64), target)
    print(count, next(read_fastx(target)).qualities)
```

- **Example output:**

```text
1 (30, 31, 32, 33)
```

[Dependencies and implementation notes](../../faq/bioinformatics.md#reads).

## 8) Read JASPAR, MEME or TRANSFAC motifs {#bio-50}

- **Purpose:** Convert motif files into native models ready for scanning.
- **Method:** Decode counts with Biopython, then construct `MotifModel` objects.
- **API:** `dnakit.bioinformatics.motifs.read_motifs(path, format='jaspar', pseudocount=0.5)`
- **Input:** A motif file; requires `bio-formats` and `bio-analysis`. This example creates a small JASPAR file.
- **Output:** Returns a tuple of `MotifModel` entries with motif names and probability matrices.
- **Example code:**

```python
from tempfile import TemporaryDirectory
from pathlib import Path
from dnakit.bioinformatics.motifs import read_motifs

with TemporaryDirectory() as folder:
    path = Path(folder) / "motif.jaspar"
    path.write_text(">MA0000.1 example\nA [ 3 0 0 ]\nC [ 0 3 0 ]\nG [ 0 0 3 ]\nT [ 0 0 0 ]\n")
    models = read_motifs(path)
    print(len(models), models[0].consensus)
```

- **Example output:**

```text
1 ACG
```

[Dependencies and implementation notes](../../faq/bioinformatics.md#optimization).

## 9) Read Sanger ABI bases, qualities and signals {#bio-55-io}

- **Purpose:** Read sequencing signals and base calls, including a quality-trimmed sequence.
- **Method:** Decode ABI with Biopython, then trim low-quality ends with the Mott algorithm.
- **API:** `dnakit.bioinformatics.cloning.read_sanger(path, error_probability_cutoff=0.05)`
- **Input:** An existing `read.ab1` file; requires `bio-formats`.
- **Output:** Returns a mapping with the trim interval, trimmed sequence, quality scores, and signals.
- **Example code:**

```python
from dnakit.bioinformatics.cloning import read_sanger

trace = read_sanger("read.ab1")
print(trace["trim_interval"], trace["trimmed_sequence"])
```

- **Example output:**

`trace` contains `sequence`, `qualities`, `trim_interval`, `trimmed_sequence`, `channels` and `peak_positions`; the example prints the trim interval and trimmed sequence. Requires an actual ABI file; this example was not executed.

[Dependencies and implementation notes](../../faq/bioinformatics.md#optimization).

## 10) Read SnapGene or GenBank constructs {#bio-57}

- **Purpose:** Read construct sequences and feature annotations.
- **Method:** Decode strictly with Biopython by default; explicit `permissive=True` uses pydna.
- **API:** `dnakit.bioinformatics.cloning.read_constructs(path, format='snapgene', permissive=False)`
- **Input:** An existing `.dna` or GenBank file; strict mode requires `bio-formats`, permissive mode requires `bio-assembly`.
- **Output:** Returns a tuple of construct records with sequences and features.
- **Example code:**

```python
from dnakit.bioinformatics.cloning import read_constructs

records = read_constructs("construct.gb", format="genbank")
for record in records:
    print(record.id, len(record.seq), len(record.features))
```

- **Example output:**

Returns a tuple of records; the loop prints each construct ID, sequence length and feature count. Requires an actual construct file; this example was not executed.

[Dependencies and implementation notes](../../faq/bioinformatics.md#optimization).

## Convenient row and result persistence {#streaming-row-artifacts}

Use these `dnakit.io` functions for ordinary result files. Existing `read_table` and `export_table` retain their explicit typed schema and resource-limit contracts.

| API | Purpose |
| --- | --- |
| `iter_rows(source, batch_size=4096)` | Stream CSV, TSV, JSONL and Parquet rows; JSON arrays are loaded in memory |
| `write_rows(rows, target, columns=None, overwrite=False, batch_size=4096)` | Consume rows once, write one path or several paths, and return the row count |
| `RowWriter(target, columns=None, overwrite=False, batch_size=4096)` | Append `.write(row)` inside a context manager, including routing rows to separate groups |
| `save_json(value, target, overwrite=False)` / `load_json(source)` | Preserve nested JSON objects and arrays; saving also accepts results with `to_dict()` |
| `save_arrays(arrays, target, compressed=True, overwrite=False)` / `load_arrays(source)` | Persist named NPZ arrays without pickle; object arrays are rejected |
| `file_digest(source, algorithm="sha256")` | Compute file checksums with bounded memory |

The target suffix selects CSV, TSV, JSON array, JSONL, Parquet or FASTA. Text formats support `.gz`. Columns are inferred from the first row or supplied explicitly; subsequent rows must have exactly those keys. Empty collections require columns. CSV/TSV write `None` as an empty cell and read values as strings; JSON/JSONL/Parquet preserve their own types. FASTA uses only `id` and `sequence`; include a tabular target to retain other fields. Continue using `dnakit.read()` to read FASTA.

```python
from dnakit.io import iter_rows, write_rows, save_json

rows = ({"id": f"seq{i}", "sequence": sequence, "label": i}
        for i, sequence in enumerate(["ACGT", "NNRY"]))
count = write_rows(rows, ["results/sequences.csv", "results/sequences.fasta"], overwrite=True)
save_json({"count": count, "parameters": {"source": "demo"}}, "results/run.json", overwrite=True)
print(list(iter_rows("results/sequences.csv")))
```

Writers create parents, reject overwriting by default, and stage files until successful completion. Each final replacement is atomic; several replacements are not a filesystem transaction. Text row writers do not retain the full collection. Parquet buffers up to `batch_size` rows and fixes types from the first batch; use `TableSchema` with `export_table` for explicit types. Parquet requires PyArrow and array persistence requires NumPy. `load_arrays` loads arrays into memory and closes the archive before returning.
