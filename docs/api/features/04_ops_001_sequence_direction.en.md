# OPS-001 sequence completion

Public processing APIs share [file and text input rules](00_processing_inputs.md).

Reverse, complement, or reverse complement the DNA sequence to obtain a new sequence after the direction is reversed.

## 1) OPS-001.1 reverse order

- **Function:** Only reverse the character order of the DNA sequence without replacing the bases, used to check or construct the reverse position sequence.
- **API:** `dnakit.ops.reverse(sequence[required])`.
- **Input:** One required `DNASequence`; no additional required parameters.
- **Sample code:**

```python
from dnakit import DNASequence
from dnakit.ops import reverse

seq = DNASequence("AATCGC")
print(reverse(seq).symbols)
```

- **Example results:**

```text
CGCTAA
```

## 2) OPS-001.2 complementary

- **Function:** Replace each DNA base with a complementary base according to base pairing rules, but maintain the original position order, and use it to construct the corresponding complementary chain.
- **API:** `dnakit.ops.complement(sequence[required])`.
- **Input:** One required `DNASequence`; no additional required parameters.
- **Sample code:**

```python
from dnakit import DNASequence
from dnakit.ops import complement

seq = DNASequence("AATCGC")
print(complement(seq).symbols)
```

- **Example results:**

```text
TTAGCG
```

## 3) OPS-001.3 Reverse complementation

Large uncompressed FASTA/FNA inputs (at least 16 MiB) with default reader settings
use file-region processing. Workers read, validate, transform, and write disjoint
regions directly; only paths and offsets cross process boundaries. Record order is
preserved and the temporary output is installed only after successful validation.
`workers=1` runs the same block algorithm serially; `workers=-1` uses available CPUs
subject to memory limits. Multiprocessing scripts need an `if __name__ == "__main__":`
guard. Irregular wrapping, compressed files, other formats, explicit `read_config`,
and output byte limits retain the general processing path and its validation rules.

- **Function:** Simultaneously reverse the sequence order and replace it with complementary bases to obtain a 5′→3′ representation of the other strand, which can be used for reverse strand search, direction unification and primer analysis.
- **API:** `dnakit.ops.reverse_complement(value, output_path=None, *, feature_policy="preserve", read_config=None, write_config=None, sequence_column=None, show_progress=True)`.
- **Input:** DNA text, a file path, or an existing single-record `DNA`, `DNASequence`, or `DNARecord`. Text returns a reverse-complement string; existing objects retain their return types while mapping features, strands, and per-base annotations. File suffixes identify `.fna/.fa/.fasta`, `.txt`, `.csv/.tsv`, and other DNAKit I/O formats, including `.gz`. Strings containing a dot or path separator are treated as paths; use `Path` explicitly for extensionless files.
- **Sample code:**

```python
from dnakit import DNA
from dnakit.ops import reverse_complement

dna = DNA("AATCGC", id="seq-1")
print(reverse_complement(dna).symbols)
```

- **Example results:**

```text
GCGATT
```

DNA strings and sequence files are also supported directly:

```python
from dnakit.ops import reverse_complement

print(reverse_complement("AATCGC"))  # GCGATT; accepts lowercase and whitespace
result = reverse_complement("hg38.fna")  # writes hg38_reverse_complement.fna
result = reverse_complement("sequences.txt", "rc.txt")
result = reverse_complement("sequences.csv", "rc.csv", sequence_column="sequence")
```

Strings containing a dot or path separator are treated as paths; use `Path` for extensionless files.
Files infer FASTA/FNA, TXT, CSV/TSV and other DNAKit I/O formats, including gzip, and stream one record
at a time. The result is a `WriteResult` with `record_count` and `target_artifact.relative_path`
(the output path relative to the current working directory).
File mode displays start/completion messages and processed record counts; set `show_progress=False`
to silence it. DNA strings return uppercase IUPAC strings; existing DNA object return types remain unchanged.

Without `read_config`, uncompressed input limits follow the actual file size, so hg38 requires no
manual size settings. Gzip input has a 20 GB decompressed reading ceiling. Pass `read_config=ReadConfig(...)` to set explicit limits. Input and output must differ; existing outputs require
`write_config=WriteConfig(overwrite=True)` to replace them.

TXT defaults to one sequence per non-empty line and writes one sequence per line without IDs.
Use `read_config=ReadConfig(text_mode="sequence")` for one wrapped sequence spanning the whole file.
CSV reads the `sequence` column by default, with `sequence_column` allowing another name, and preserves
the input column names, order, labels and other extra columns. Record annotations follow the reverse-complement coordinate mapping.
