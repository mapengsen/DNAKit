# Optional sequence backends: Biopython, scikit-bio, pyfastx, screed and py2bit

This group of functions in `dnakit.sequence` is a set of adapters over optional packages. EMBL and GenBank records are parsed, extracted and converted with [Biopython](https://github.com/biopython/biopython)'s `Bio.SeqIO`. DNA translation, k-mer frequencies and degenerate-base expansion use [scikit-bio](https://github.com/scikit-bio/scikit-bio). FASTQ quality encodings are inferred with [pyfastx](https://github.com/lmdu/pyfastx). Record databases are written and queried with [screed](https://github.com/dib-lab/screed). Region reads and summaries of 2bit files use [py2bit](https://github.com/deeptools/py2bit). Every function works only on an explicit file or sequence. None of them accesses the network, and none computes primers, alignments or restriction digests.

The dependencies are optional and are imported only when the matching function is called. If a package is missing, the function raises an error that includes the install command.

## Example data

The examples below import the functions directly from `dnakit.sequence` and use these two files. The EMBL file has 45 bases and two CDS features. The FASTQ file has two short reads. Both are only for demonstrating the interface. The 2bit file is binary, so the example code writes it itself.

**Example file:** `records.embl`

```text
ID   AB000001; SV 1; linear; genomic DNA; STD; UNC; 45 BP.
XX
DE   test record for dnakit
XX
FH   Key             Location/Qualifiers
FH
FT   CDS             1..18
FT                   /gene="geneA"
FT   CDS             complement(21..32)
FT                   /gene="geneB"
XX
SQ   Sequence 45 BP;
     atggccaaat ttgggcccta aggatccatg aaattttaac cggta                        45
//
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

## 1) `dnakit.sequence.read_sequence_records()`

- **Purpose:** Parse EMBL or GenBank records and return each record's sequence, length and feature counts.
- **Method:** Read records one at a time with `Bio.SeqIO.parse`. If the record count or base count exceeds `max_records` or `max_bases`, the function raises an error instead of truncating the result.
- **API:** `dnakit.sequence.read_sequence_records(source, *, format="embl", max_records=1000, max_bases=10000000)`.
- **Input:** An EMBL or GenBank file path. `format` is `"embl"` or `"genbank"`.
- **Output:** A dictionary with `record_count`, `total_bases` and `records`. Each record has `id`, `name`, `description`, `length`, `feature_counts` and `sequence`.
- **Example code:**

```python
from dnakit.sequence import read_sequence_records

result = read_sequence_records("records.embl", format="embl")
record = result["records"][0]
print(record["id"], record["length"], record["feature_counts"])
```

- **Example output:**

```text
AB000001.1 45 {'CDS': 2}
```

## 2) `dnakit.sequence.extract_feature_sequences()`

- **Purpose:** Extract sequences of one feature type, such as CDS or gene.
- **Method:** Call `extract()` from `Bio.SeqFeature`. Minus-strand features are reverse-complemented, and `join()` locations are concatenated from their parts. The output `start` and `end` cover the whole feature as a 1-based inclusive interval.
- **API:** `dnakit.sequence.extract_feature_sequences(source, *, format="genbank", feature_type="CDS", max_records=1000, max_features=10000, max_bases=10000000)`.
- **Input:** An EMBL or GenBank file path. `feature_type` is the feature type name and cannot contain whitespace.
- **Output:** `feature_count` and `features`. Each feature has `record_id`, `feature_index`, `gene`, `locus_tag`, `start`, `end`, `strand`, `parts`, `length` and `sequence`.
- **Example code:**

```python
from dnakit.sequence import extract_feature_sequences

result = extract_feature_sequences("records.embl", format="embl", feature_type="CDS")
for feature in result["features"]:
    print(feature["gene"], feature["start"], feature["end"], feature["sequence"])
```

- **Example output:**

```text
geneA 1 18 ATGGCCAAATTTGGGCCC
geneB 21 32 TTCATGGATCCT
```

## 3) `dnakit.sequence.convert_sequence_file()`

- **Purpose:** Convert EMBL or GenBank records to FASTA, EMBL or GenBank.
- **Method:** Write records with `Bio.SeqIO.write`. The target file must not exist. If the record count exceeds `max_records`, the partly written file is deleted. FASTA output keeps only the header and the sequence.
- **API:** `dnakit.sequence.convert_sequence_file(source, target, *, input_format="genbank", output_format="fasta", max_records=1000)`.
- **Input:** The source file path, a new target path, the input format (`"embl"` or `"genbank"`) and the output format (`"fasta"`, `"embl"` or `"genbank"`).
- **Output:** `record_count`, `source`, `target`, `input_format` and `output_format`.
- **Example code:**

```python
import os

from dnakit.sequence import convert_sequence_file

if os.path.exists("records.fasta"):
    os.remove("records.fasta")  # conversion never overwrites a file
result = convert_sequence_file(
    "records.embl", "records.fasta", input_format="embl", output_format="fasta"
)
with open("records.fasta", encoding="utf-8") as handle:
    print(result["record_count"], handle.readline().strip())
```

- **Example output:**

```text
1 >AB000001.1 test record for dnakit
```

## 4) `dnakit.sequence.translate_dna_with_genetic_code()`

- **Purpose:** Translate DNA with a chosen NCBI genetic code.
- **Method:** Call scikit-bio's `DNA.translate`, which transcribes to RNA and then translates. `reading_frame` is 1, 2 or 3 for the forward strand, or -1, -2 or -3 for the reverse complement. `start` and `stop` are `"ignore"`, `"require"` or `"optional"`. With `"require"`, the function raises an error when the required codon is absent. Sequences with gaps are rejected.
- **API:** `dnakit.sequence.translate_dna_with_genetic_code(sequence, *, genetic_code=1, reading_frame=1, start="ignore", stop="ignore", max_length=1000000)`.
- **Input:** A DNA string that may contain IUPAC ambiguity codes. `genetic_code` is an NCBI table ID supported by scikit-bio, such as 1, 2 or 11.
- **Output:** `protein`, `protein_length`, `genetic_code` (with `table_id` and `name`), `reading_frame`, `start` and `stop`.
- **Example code:**

```python
from dnakit.sequence import translate_dna_with_genetic_code

result = translate_dna_with_genetic_code("ATGTGAATAAGA", genetic_code=2)
print(result["protein"], result["genetic_code"]["name"])
```

- **Example output:**

```text
MWM* Vertebrate Mitochondrial
```

## 5) `dnakit.sequence.calculate_kmer_frequencies()`

- **Purpose:** Count k-mers as raw counts or relative frequencies.
- **Method:** Call scikit-bio's `DNA.kmer_frequencies`. With `overlap=False`, k-mers do not overlap. With `relative=True`, each value is its count divided by the total number of k-mers counted. k-mers that contain N or other ambiguity codes are counted as written.
- **API:** `dnakit.sequence.calculate_kmer_frequencies(sequence, *, k=3, overlap=True, relative=False, max_length=1000000)`.
- **Input:** A DNA string. `k` is an integer from 1 to 32.
- **Output:** `kmers` (a dictionary sorted by k-mer), `distinct_kmers`, `k`, `overlap`, `relative` and `sequence_length`.
- **Example code:**

```python
from dnakit.sequence import calculate_kmer_frequencies

result = calculate_kmer_frequencies("ACACATTTATTA", k=3, overlap=False, relative=True)
print(result["kmers"])
```

- **Example output:**

```text
{'ACA': 0.25, 'CAT': 0.25, 'TTA': 0.5}
```

## 6) `dnakit.sequence.expand_degenerate_sequences()`

- **Purpose:** List every definite DNA sequence represented by a degenerate sequence.
- **Method:** Call scikit-bio's `DNA.expand_degenerates`. Before expansion, the function computes the result count, which is the product of the number of allowed bases at each degenerate position. If that count exceeds `max_sequences`, the function raises an error and generates nothing. Results are in lexicographic order.
- **API:** `dnakit.sequence.expand_degenerate_sequences(sequence, *, max_sequences=1024, max_length=1000000)`.
- **Input:** An IUPAC DNA string, such as a sequence containing R, Y or N.
- **Output:** `sequences`, `expansion_count`, `degenerate_positions` and `sequence_length`.
- **Example code:**

```python
from dnakit.sequence import expand_degenerate_sequences

print(expand_degenerate_sequences("TRG")["sequences"])
```

- **Example output:**

```text
['TAG', 'TGG']
```

## 7) `dnakit.sequence.detect_fastq_quality_encoding()`

- **Purpose:** Infer the FASTQ quality offset (Phred+33 or Phred+64) and the quality score range.
- **Method:** Call pyfastx's `Fastq` without building an index. The function scans the whole file and reads the offset and character range that pyfastx reports. When pyfastx returns 0, both 33 and 64 fit the observed characters. In that case `phred_offset` is `None` and no score range is given. The file size is checked against `max_bytes` before the scan.
- **API:** `dnakit.sequence.detect_fastq_quality_encoding(source, *, max_bytes=1000000000)`.
- **Input:** A FASTQ file path. Gzip-compressed files are accepted.
- **Output:** `phred_offset`, `possible_phred_offsets`, `min_quality_ascii`, `max_quality_ascii`, `min_quality_score`, `max_quality_score`, `min_length` and `max_length`.
- **Example code:**

```python
from dnakit.sequence import detect_fastq_quality_encoding

result = detect_fastq_quality_encoding("reads.fq")
print(result["phred_offset"], result["min_quality_score"], result["max_quality_score"])
```

- **Example output:**

```text
33 0 40
```

## 8) `dnakit.sequence.build_screed_database()`

- **Purpose:** Write a FASTA or FASTQ file to a new screed SQLite record database.
- **Method:** Call `screed.create_db`, with field definitions taken from screed's FASTA or FASTQ reader. If the path does not end in `_screed`, that suffix is added. If the target already exists, the function raises an error. If the record count exceeds `max_records`, the unfinished database is deleted.
- **API:** `dnakit.sequence.build_screed_database(source, *, database, max_records=1000)`.
- **Input:** A FASTA or FASTQ file path and a database path that does not exist yet.
- **Output:** `database` (the path actually written), `format` (`"fasta"` or `"fastq"`) and `record_count`.
- **Example code:**

```python
import os

from dnakit.sequence import build_screed_database

if os.path.exists("reads_screed"):
    os.remove("reads_screed")  # building never overwrites a database
result = build_screed_database("reads.fq", database="reads")
print(result["format"], result["record_count"])
```

- **Example output:**

```text
fastq 2
```

## 9) `dnakit.sequence.fetch_screed_records()`

- **Purpose:** Read named records from a screed database. Without names, read all records.
- **Method:** Look up keys with `screed.ScreedDB`. Sequences are converted to uppercase. FASTQ records also return `quality`. Missing names raise an error that lists up to 10 of them.
- **API:** `dnakit.sequence.fetch_screed_records(database, *, names=None, max_records=1000)`.
- **Input:** The database path returned by `build_screed_database()`, or the same path without the `_screed` suffix. `names` is a list of record names.
- **Output:** `record_count` and `records`. Each record has `name`, `description`, `length` and `sequence`. FASTQ records also have `quality`.
- **Example code:**

```python
from dnakit.sequence import fetch_screed_records

record = fetch_screed_records("reads", names=["q1"])["records"][0]
print(record["name"], record["sequence"], record["quality"])
```

- **Example output:**

```text
q1 ACGTAC IIII#!
```

## 10) `dnakit.sequence.summarize_twobit()`

- **Purpose:** Summarize each sequence in a 2bit file: length, hard-masked (N) bases, ACGT bases and GC content.
- **Method:** Use py2bit's `chroms()`, `bases()` and `hardMaskedBlocks()`. Soft-masked (lowercase) bases are ignored. The GC percentage uses ACGT bases as its denominator, and N bases are not counted.
- **API:** `dnakit.sequence.summarize_twobit(source, *, max_chromosomes=10000)`.
- **Input:** A 2bit file path. The example below first writes a small file, `toy.2bit`, that contains only ACGT bases.
- **Output:** `chromosome_count`, `total_length`, `total_hard_masked` and `chromosomes`. Each sequence has `name`, `length`, `acgt_bases`, `hard_masked_bases` and `gc_percent_of_acgt`.
- **Example code:**

```python
import struct

from dnakit.sequence import summarize_twobit

# Minimal little-endian .2bit writer for this example: ACGT only, no masks.
records = [("chr1", "ACGTACGTGGCC"), ("chr2", "TTTTAAAACCCCGGGG")]
codes = {"T": 0, "C": 1, "A": 2, "G": 3}
index, body = b"", b""
offset = 16 + sum(1 + len(name) + 4 for name, _ in records)
for name, sequence in records:
    packed = bytearray((len(sequence) + 3) // 4)
    for position, base in enumerate(sequence):
        packed[position // 4] |= codes[base] << (6 - 2 * (position % 4))
    block = struct.pack("<IIII", len(sequence), 0, 0, 0) + bytes(packed)
    index += bytes([len(name)]) + name.encode() + struct.pack("<I", offset)
    body += block
    offset += len(block)
with open("toy.2bit", "wb") as handle:
    handle.write(struct.pack("<IIII", 0x1A412743, 0, len(records), 0) + index + body)

result = summarize_twobit("toy.2bit")
print([(c["name"], c["length"], c["gc_percent_of_acgt"]) for c in result["chromosomes"]])
```

- **Example output:**

```text
[('chr1', 12, 66.66666666666667), ('chr2', 16, 50.0)]
```

## 11) `dnakit.sequence.fetch_twobit_region()`

- **Purpose:** Read one interval of one sequence from a 2bit file.
- **Method:** Call py2bit's `sequence()`. Coordinates are 0-based and half-open. If `end` is beyond the sequence length, the function raises an error instead of silently truncating the result as py2bit does. The output is uppercase.
- **API:** `dnakit.sequence.fetch_twobit_region(source, chromosome, start, end, *, max_length=1000000)`.
- **Input:** A 2bit file path, a sequence name, a 0-based `start` and an exclusive `end`. The example reads the `toy.2bit` file written in the previous section.
- **Output:** `chromosome`, `start`, `end`, `length`, `hard_masked_bases` and `sequence`.
- **Example code:**

```python
from dnakit.sequence import fetch_twobit_region

print(fetch_twobit_region("toy.2bit", "chr1", 2, 8)["sequence"])
```

- **Example output:**

```text
GTACGT
```

## Installing the optional dependencies

```bash
python -m pip install "dnakit[bio-formats]"     # biopython: EMBL/GenBank reading, extraction and conversion
python -m pip install "dnakit[bio-skibio]"      # scikit-bio: translation, k-mers and degenerate-base expansion
python -m pip install "dnakit[bio-pyfastx]"     # pyfastx: FASTQ quality encoding
python -m pip install "dnakit[bio-screed]"      # screed: SQLite record database
python -m pip install "dnakit[bio-py2bit]"      # py2bit: 2bit summaries and region reads
```

Project pages: [Biopython](https://github.com/biopython/biopython), [scikit-bio](https://github.com/scikit-bio/scikit-bio), [pyfastx](https://github.com/lmdu/pyfastx), [screed](https://github.com/dib-lab/screed), [py2bit](https://github.com/deeptools/py2bit).

## How these functions are verified

The test file is `others/tests/unit/test_sequence_backends.py`, with 32 tests. They check the following:

- EMBL/GenBank reading, feature extraction and format conversion: gene sequences, minus-strand reverse complements and `join()` locations match hand calculations. The input files are generated by Biopython inside the tests.
- scikit-bio wrappers: translation with different genetic codes and reading frames, start and stop policies, non-overlapping and relative k-mer counts, and the count and content of degenerate expansions all match hand calculations. Overlapping k-mer counts also match a direct scikit-bio call entry by entry.
- pyfastx, screed and py2bit: FASTQ quality offsets and score ranges, screed FASTA and FASTQ round trips, 2bit GC and hard-mask counts, and 0-based half-open region reads all match hand calculations.
- Exceeding `max_records`, `max_bases`, `max_features`, `max_length`, `max_sequences` or `max_bytes` raises an error. A failed conversion or database build leaves no unfinished output file.
- A missing optional package raises an install hint. Gap characters, unknown genetic-code IDs, files that are not 2bit, inputs that are not FASTQ, and output paths that already exist are all rejected.

## Current scope

This batch wraps only DNA and gene-feature operations that dnakit does not already provide in the same form. The following operations are not wrapped:

- pyfaidx region fetching is already provided by `dnakit.integrations.fetch_pyfaidx_region()`.
- FASTA/FASTQ summary statistics from pyfastx are already provided by `dnakit.integrations.summarize_pyfastx()`, and streaming summaries from dnaio by `dnakit.integrations.summarize_dnaio_fastx()`. Writing, pairing and reverse-complementing dnaio records are provided natively by `dnakit.bioinformatics.reads` and `dnakit.fasta_fastq`.
- FASTA/FASTQ parsing with Biopython `Bio.SeqIO` overlaps with the native `dnakit.bioinformatics.reads.read_fastx()`, so it is not wrapped. GenBank parsing overlaps with the native `dnakit.io.genbank` module. It is kept here because feature extraction and format conversion share the Biopython path with EMBL parsing.
- Reverse complement and GC content are provided natively by `complement_dna`, `reverse_complement_dna` and `calculate_sequence_metrics`, so scikit-bio's equivalent methods are not wrapped.
- The native `translate_dna_to_protein()` supports only the standard genetic code. The translation function in this batch uses scikit-bio, which supports several genetic codes. Biopython's `Seq.translate` would duplicate it, so it is not wrapped.
- Melting temperature from Biopython `MeltingTemp` overlaps with the native `dnakit.thermodynamics.melting_temperature()`, so it is not wrapped.
- Hamming distance is provided natively by `dnakit.similarity.hamming_distance()`, so scikit-bio's `hamming` is not wrapped.
- Restriction digestion is provided natively by `dnakit.patterns.scan_restriction_sites()` and `dnakit.molbio.digest_restriction()`, so Biopython's `Bio.Restriction` is not wrapped.
- GC skew is already provided by the native descriptors and is not wrapped.
- py2bit is read-only and has no writing interface, so this batch has no 2bit writing function.
- pyfastx random access needs a `.fxi` index written next to the input file. This batch does not use that feature.
