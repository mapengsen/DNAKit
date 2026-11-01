# Biopython sequence formats and Seq helpers: QUAL, IMGT, GFA, nib, PIR, SeqXML, SFF, Xdna, Entrez XML and BGZF offsets

This group of functions in `dnakit.sequence.biopython_p18_functions` wraps [Biopython](https://biopython.org/) (`Bio.SeqIO`, `Bio.Seq`, `Bio.SeqFeature`, `Bio.Entrez`, `Bio.bgzf`, `Bio.UniProt.GOA`, `Bio.Emboss`, `Bio.Align.substitution_matrices`) for the operations that the earlier Biopython wrappers did not reach: reading and writing the less common sequence file formats, exact-match operations on sequence strings, INSDC location strings, EMBOSS `primersearch` input, GO-annotation filtering, substitution-matrix text, saved NCBI Entrez XML files and BGZF virtual offsets.

Readers read local files. Writers never overwrite: the output path must be a new file. Sequence-string functions keep the case and the characters of the input, as `Bio.Seq` does. No function uses the network; `parse_entrez_xml()` stops instead of downloading a missing DTD. Biopython is an optional dependency (`pip install "dnakit[bio-formats]"`) imported only when a function is called; a missing package raises an error that includes the install command. The examples import from the module, `dnakit.sequence.biopython_p18_functions`.

## 1) `dnakit.sequence.biopython_p18_functions.read_sequence_file()`

- **Purpose:** Read the records of a less common `Bio.SeqIO` format and return them as plain dictionaries.
- **Method:** Call `Bio.SeqIO.parse()` for the chosen format: `qual` (PHRED qualities only; the sequence is undefined and returned as null), `imgt`, `genbank-cds` and `embl-cds` (one record per CDS feature, holding the translation), `gck`, `gfa1`, `gfa2` (assembly-graph segments), `ig`, `nib`, `pir`, `seqxml`, `tab`, `xdna`, `abi` and `abi-trim` (Biopython's Mott-style trimming with error cutoff 0.05). The raw ABIF block of ABI files is left out unless `include_raw_trace=True`. The call stops when the file holds more than `max_records` records.
- **API:** `dnakit.sequence.biopython_p18_functions.read_sequence_file(path, *, file_format='tab', max_records=10000, include_sequence=True, include_raw_trace=False)`.
- **Input:** A local file path, `file_format`, optional `max_records` (default 10000), `include_sequence` and `include_raw_trace`.
- **Output:** `record_count` and `records`; each record has `id`, `name`, `description`, `length`, `sequence`, `annotations`, `letter_annotations` (for QUAL the `phred_quality` list), `feature_count`, `features` (type, location, qualifiers) and `dbxrefs`.
- **Example code:**

```python
from pathlib import Path

from dnakit.sequence.biopython_p18_functions import read_sequence_file

Path("records.tab").write_text("seq1\tACGTACGTAC\nseq2\tGGTTAACC\n", encoding="utf-8")
result = read_sequence_file("records.tab", file_format="tab")
print(result["record_count"], [(r["id"], r["length"]) for r in result["records"]])
print(result["records"][0]["sequence"])
```

- **Example output:**

```text
2 [('seq1', 10), ('seq2', 8)]
ACGTACGTAC
```

## 2) `dnakit.sequence.biopython_p18_functions.write_sequence_file()`

- **Purpose:** Write sequence records to a new file in an IMGT, nib, PHD, PIR, SeqXML, tab or Xdna format.
- **Method:** Build `SeqRecord` objects from the mappings and call `Bio.SeqIO.write()`. The `nib` and `xdna` formats hold exactly one record. SFF files carry flow data, so they are written by `convert_sequence_format()` from an SFF source instead. The output path must not exist; a failed write removes the partial file.
- **API:** `dnakit.sequence.biopython_p18_functions.write_sequence_file(records, output_path, *, file_format='tab')`.
- **Input:** A list of record mappings with `id` and `sequence` and, optionally, `name`, `description`, `molecule_type` (`DNA` or `RNA`), `topology` (`linear` or `circular`) and `qualities` (PHRED integers, required for `phd`); a new output path; `file_format`.
- **Output:** `output_path`, `file_format`, `record_count` and `bytes_written`.
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.sequence.biopython_p18_functions import write_sequence_file

out = Path(tempfile.mkdtemp()) / "demo.pir"
records = [
    {"id": "a", "sequence": "ACGTACGT", "description": "first"},
    {"id": "b", "sequence": "GGTTAACC"},
]
result = write_sequence_file(records, out, file_format="pir")
print(result["record_count"], result["bytes_written"])
print(out.read_text(encoding="utf-8"))
```

- **Example output:**

```text
2 44
>D1;a
a - first
ACGTACGT*
>D1;b
b
GGTTAACC*
```

## 3) `dnakit.sequence.biopython_p18_functions.convert_sequence_format()`

- **Purpose:** Convert a sequence file between `Bio.SeqIO` formats, including SFF output.
- **Method:** Parse the source with `Bio.SeqIO.parse()` and write with `Bio.SeqIO.write()`. Source formats are those of `read_sequence_file()` plus embl, fasta, fastq, genbank, sff and sff-trim; targets are embl, fasta, fastq, genbank, imgt, nib, phd, pir, qual, seqxml, sff, tab and xdna. For an SFF target the writer is `SffWriter`, which copies the Roche XML manifest of an SFF source. Information a target cannot hold (features in FASTA, qualities in tab) is dropped, as Biopython does.
- **API:** `dnakit.sequence.biopython_p18_functions.convert_sequence_format(source, output_path, *, source_format='genbank', target_format='tab', max_records=10000)`.
- **Input:** A source path, a new output path, `source_format`, `target_format` and `max_records`.
- **Output:** `record_count`, `source_format`, `target_format`, `roche_manifest_copied` and the paths.
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.sequence.biopython_p18_functions import convert_sequence_format

Path("records.tab").write_text("seq1\tACGTACGTAC\nseq2\tGGTTAACC\n", encoding="utf-8")
out = Path(tempfile.mkdtemp()) / "records.fasta"
result = convert_sequence_format("records.tab", out, source_format="tab", target_format="fasta")
print(result["record_count"], result["target_format"])
print(out.read_text(encoding="utf-8"))
```

- **Example output:**

```text
2 fasta
>seq1
ACGTACGTAC
>seq2
GGTTAACC
```

## 4) `dnakit.sequence.biopython_p18_functions.read_sff_manifest()`

- **Purpose:** Return the Roche XML manifest stored in an SFF file.
- **Method:** Call `Bio.SeqIO.SffIO.ReadRocheXmlManifest()` on the binary file. A file without a manifest gives `has_manifest = False`; a file that is not SFF raises an error.
- **API:** `dnakit.sequence.biopython_p18_functions.read_sff_manifest(path)`.
- **Input:** A local SFF file path.
- **Output:** `has_manifest`, `xml` and `xml_length`.
- **Example code:**

```python
from Bio.Seq import Seq
from Bio.SeqIO import SffIO
from Bio.SeqRecord import SeqRecord

from dnakit.sequence.biopython_p18_functions import read_sff_manifest

record = SeqRecord(Seq("ACGT"), id="read1")
record.letter_annotations["phred_quality"] = [30, 31, 32, 33]
record.annotations.update(
    flow_values=(0, 100, 100, 100, 100, 0, 0, 0), flow_index=(2, 1, 1, 1), flow_chars="TACGTACG",
    flow_key="", clip_qual_left=1, clip_qual_right=4, clip_adapter_left=0, clip_adapter_right=0,
)
with open("demo.sff", "wb") as handle:
    SffIO.SffWriter(handle, xml="<manifest><run/></manifest>").write_file([record])
result = read_sff_manifest("demo.sff")
print(result["has_manifest"], result["xml"], result["xml_length"])
```

- **Example output:**

```text
True <manifest><run/></manifest> 27
```

## 5) `dnakit.sequence.biopython_p18_functions.count_subsequence()`

- **Purpose:** Count the occurrences of an exact subsequence, non-overlapping or overlapping.
- **Method:** `overlap=False` calls `Seq.count()` (non-overlapping, scanning left to right); `overlap=True` calls `Seq.count_overlap()`. The match is case sensitive. `start` and `end` are optional 0-based slice limits.
- **API:** `dnakit.sequence.biopython_p18_functions.count_subsequence(sequence, subsequence, *, overlap=False, start=None, end=None)`.
- **Input:** `sequence`, `subsequence`, `overlap`, `start`, `end`.
- **Output:** `count`, `overlap`, `sequence_length` and `subsequence`.
- **Example code:**

```python
from dnakit.sequence.biopython_p18_functions import count_subsequence

sequence = "AAAAGAAAAGAAAA"
print(count_subsequence(sequence, "AA")["count"])
print(count_subsequence(sequence, "AA", overlap=True)["count"])
print(count_subsequence(sequence, "AA", overlap=True, start=2, end=8)["count"])
```

- **Example output:**

```text
6
9
3
```

## 6) `dnakit.sequence.biopython_p18_functions.search_subsequences()`

- **Purpose:** Find every position of several exact subsequences.
- **Method:** Call `Seq.search()`, which reports each position of each subsequence; matches may overlap each other. Positions are 0-based and in increasing order; the match is case sensitive.
- **API:** `dnakit.sequence.biopython_p18_functions.search_subsequences(sequence, subsequences)`.
- **Input:** `sequence` and a list of up to 1000 `subsequences`.
- **Output:** `hit_count` and `hits`, each with `position` and `subsequence`.
- **Example code:**

```python
from dnakit.sequence.biopython_p18_functions import search_subsequences

result = search_subsequences("AAAAGAAAAGAAAA", ["AAG", "GA"])
print(result["hit_count"])
for hit in result["hits"]:
    print(hit["position"], hit["subsequence"])
```

- **Example output:**

```text
4
2 AAG
4 GA
7 AAG
9 GA
```

## 7) `dnakit.sequence.biopython_p18_functions.split_sequence()`

- **Purpose:** Split a sequence at a separator.
- **Method:** Call `Seq.split()`, or `Seq.rsplit()` when `from_right=True`. Without a separator the split is on whitespace. `maxsplit` limits the number of splits (-1 means no limit).
- **API:** `dnakit.sequence.biopython_p18_functions.split_sequence(sequence, separator=None, *, maxsplit=-1, from_right=False)`.
- **Input:** `sequence`, optional `separator`, `maxsplit`, `from_right`.
- **Output:** `piece_count` and `pieces`.
- **Example code:**

```python
from dnakit.sequence.biopython_p18_functions import split_sequence

sequence = "AAAAGAAAAGAAAA"
print(split_sequence(sequence, "G")["pieces"])
print(split_sequence(sequence, "G", maxsplit=1)["pieces"])
print(split_sequence(sequence, "G", maxsplit=1, from_right=True)["pieces"])
```

- **Example output:**

```text
['AAAA', 'AAAA', 'AAAA']
['AAAA', 'AAAAGAAAA']
['AAAAGAAAA', 'AAAA']
```

## 8) `dnakit.sequence.biopython_p18_functions.trim_sequence()`

- **Purpose:** Trim the ends of a sequence the way `Bio.Seq.Seq` does.
- **Method:** `mode` is `strip` (both ends), `lstrip` or `rstrip`, which remove any character of `pattern` (whitespace when it is null), or `removeprefix` / `removesuffix`, which remove the exact `pattern` string once. Case is kept.
- **API:** `dnakit.sequence.biopython_p18_functions.trim_sequence(sequence, *, mode='strip', pattern=None)`.
- **Input:** `sequence`, `mode` and `pattern`.
- **Output:** `trimmed`, `trimmed_length`, `original_length`, `mode` and `pattern`.
- **Example code:**

```python
from dnakit.sequence.biopython_p18_functions import trim_sequence

for mode, pattern in [("strip", "A"), ("lstrip", "A"), ("removeprefix", "AAAA"), ("removesuffix", "AAAA")]:
    print(mode, trim_sequence("AAAAGAAAAGAAAA", mode=mode, pattern=pattern)["trimmed"])
```

- **Example output:**

```text
strip GAAAAG
lstrip GAAAAGAAAA
removeprefix GAAAAGAAAA
removesuffix AAAAGAAAAG
```

## 9) `dnakit.sequence.biopython_p18_functions.parse_feature_location()`

- **Purpose:** Parse an INSDC feature-location string such as `complement(join(1..10,20..30))`.
- **Method:** Call `Bio.SeqFeature.Location.fromstring()`. The overall start is 0-based and the end exclusive; the input string is 1-based. Fuzzy positions (`<1`, `>200`, `(5.9)`, `12^13`) are reported with their Biopython position type.
- **API:** `dnakit.sequence.biopython_p18_functions.parse_feature_location(location)`.
- **Input:** A location string.
- **Output:** `location_type`, `start` and `end` (type, text, value), `strand`, `operator` (`join`, `order` or `single`), `ref`, `part_count`, `parts` and the normalized text.
- **Example code:**

```python
from dnakit.sequence.biopython_p18_functions import parse_feature_location

joined = parse_feature_location("complement(join(1..10,20..30))")
print(joined["start"]["value"], joined["end"]["value"], joined["strand"], joined["operator"])
print([(p["start"]["value"], p["end"]["value"]) for p in joined["parts"]])
fuzzy = parse_feature_location("<1..>200")
print(fuzzy["start"]["type"], fuzzy["end"]["type"], fuzzy["normalized"])
print(parse_feature_location("AB000123.1:5..10")["ref"])
```

- **Example output:**

```text
0 30 -1 join
[(19, 30), (0, 10)]
BeforePosition AfterPosition [<0:>200](+)
AB000123.1
```

## 10) `dnakit.sequence.biopython_p18_functions.format_primersearch_input()`

- **Purpose:** Create the input file of EMBOSS `primersearch` from named forward and reverse primers.
- **Method:** Add each pair with `Bio.Emboss.PrimerSearch.InputRecord.add_primer_set()` and return the text, one line per set: name, forward primer, reverse primer. Primers must be IUPAC DNA. When `output_path` is given the text is also written to that new file.
- **API:** `dnakit.sequence.biopython_p18_functions.format_primersearch_input(primer_sets, *, output_path=None)`.
- **Input:** A mapping of set name to `[forward, reverse]` and an optional new `output_path`.
- **Output:** `primer_set_count`, `text` and `output_path`.
- **Example code:**

```python
from dnakit.sequence.biopython_p18_functions import format_primersearch_input

result = format_primersearch_input({"p1": ["ACGTAC", "GGTACC"], "p2": ["AAAACG", "TTTTGC"]})
print(result["primer_set_count"])
print(result["text"])
```

- **Example output:**

```text
2
p1 ACGTAC GGTACC
p2 AAAACG TTTTGC
```

## 11) `dnakit.sequence.biopython_p18_functions.filter_go_annotations()`

- **Purpose:** Keep the records of a GAF annotation file whose fields match given values.
- **Method:** Read the file with `Bio.UniProt.GOA.gafiterator()` and keep a record when any one of the named fields holds any accepted value (`Bio.UniProt.GOA.record_has()`). List-valued fields such as `DB:Reference` or `Taxon_ID` match when any element does. Field names are those of Biopython's GAF records: `DB`, `DB_Object_ID`, `DB_Object_Symbol`, `Qualifier`, `GO_ID`, `DB:Reference`, `Evidence`, `With`, `Aspect`, `DB_Object_Name`, `Synonym`, `DB_Object_Type`, `Taxon_ID`, `Date`, `Assigned_By` and so on.
- **API:** `dnakit.sequence.biopython_p18_functions.filter_go_annotations(path, field_values, *, max_records=10000)`.
- **Input:** A GAF file path, a mapping of field name to a list of accepted values, and `max_records`.
- **Output:** `records_read`, `kept_count` and the kept `records`.
- **Example code:**

```python
from pathlib import Path

from dnakit.sequence.biopython_p18_functions import filter_go_annotations

row1 = ["UniProtKB", "P12345", "GENE1", "", "GO:0008150", "PMID:1", "IDA", "", "P", "gene one",
        "", "protein", "taxon:9606", "20200101", "UniProt", "", ""]
row2 = ["UniProtKB", "Q99999", "GENE2", "NOT", "GO:0003674", "PMID:2|PMID:3", "IEA", "", "F", "gene two",
        "", "protein", "taxon:10090", "20200102", "UniProt", "", ""]
Path("demo.gaf").write_text("!gaf-version: 2.1\n" + "\n".join("\t".join(r) for r in (row1, row2)) + "\n", encoding="utf-8")
result = filter_go_annotations("demo.gaf", {"DB:Reference": ["PMID:3"]})
print(result["records_read"], result["kept_count"], result["records"][0]["DB_Object_Symbol"])
print(filter_go_annotations("demo.gaf", {"Aspect": ["P"]})["records"][0]["GO_ID"])
```

- **Example output:**

```text
2 1 GENE2
GO:0008150
```

## 12) `dnakit.sequence.biopython_p18_functions.format_substitution_matrix()`

- **Purpose:** Return a bundled substitution matrix, for example `NUC.4.4`, as formatted text.
- **Method:** Load the matrix with `Bio.Align.substitution_matrices.load()` and format it with `Array.format()`. `fmt` is a printf-style number format (`%i`, `%.1f`, ...). The text starts with the matrix header comments, then the alphabet and one row per letter.
- **API:** `dnakit.sequence.biopython_p18_functions.format_substitution_matrix(name='NUC.4.4', *, fmt='%i')`.
- **Input:** A bundled matrix `name` and `fmt`.
- **Output:** `name`, `alphabet`, `fmt` and `text`.
- **Example code:**

```python
from dnakit.sequence.biopython_p18_functions import format_substitution_matrix

result = format_substitution_matrix("NUC.4.4", fmt="%i")
print(result["alphabet"])
print("\n".join(result["text"].splitlines()[-16:-12]))
```

- **Example output:**

```text
ATGCSWRYKMBVHDN
   A  T  G  C  S  W  R  Y  K  M  B  V  H  D  N
A  5 -4 -4 -4 -4  1  1 -4 -4  1 -4 -1 -1 -1 -2
T -4  5 -4 -4 -4  1 -4  1  1 -4 -1 -4 -1 -1 -2
G -4 -4  5 -4  1 -4  1 -4  1 -4 -1 -1 -4 -1 -2
```

## 13) `dnakit.sequence.biopython_p18_functions.parse_entrez_xml()`

- **Purpose:** Parse a saved NCBI Entrez (E-utilities) XML result file into plain Python data.
- **Method:** `multiple=False` calls `Bio.Entrez.read()` (one result such as ESearch, EInfo or ELink); `multiple=True` calls `Bio.Entrez.parse()` (a list of records such as EFetch or ESummary sets). Biopython needs the DTD named in the file. This function does not download it: it reads the DOCTYPE first and stops with an error when the DTD is neither bundled with Biopython nor in the user DTD directory.
- **API:** `dnakit.sequence.biopython_p18_functions.parse_entrez_xml(path, *, multiple=False, max_records=1000)`.
- **Input:** A local XML file path, `multiple` and `max_records`.
- **Output:** `dtd`, `multiple`, `record_count` and `data` (nested lists and dictionaries).
- **Example code:**

```python
from pathlib import Path

from dnakit.sequence.biopython_p18_functions import parse_entrez_xml

Path("esearch.xml").write_text(
    '<?xml version="1.0" encoding="UTF-8" ?>\n'
    '<!DOCTYPE eSearchResult PUBLIC "-//NLM//DTD esearch 20060628//EN" '
    '"https://eutils.ncbi.nlm.nih.gov/eutils/dtd/20060628/esearch.dtd">\n'
    "<eSearchResult><Count>2</Count><RetMax>2</RetMax><RetStart>0</RetStart>"
    "<IdList><Id>111</Id><Id>222</Id></IdList><TranslationSet></TranslationSet>"
    "<QueryTranslation>BRCA1[Gene]</QueryTranslation></eSearchResult>\n",
    encoding="utf-8",
)
result = parse_entrez_xml("esearch.xml")
print(result["dtd"], result["data"]["Count"], result["data"]["IdList"], result["data"]["QueryTranslation"])
```

- **Example output:**

```text
esearch.dtd 2 ['111', '222'] BRCA1[Gene]
```

## 14) `dnakit.sequence.biopython_p18_functions.bgzf_virtual_offset()`

- **Purpose:** Combine or split a BGZF virtual file offset.
- **Method:** Give `block_start` (compressed offset of the BGZF block, below 2**48) and `within_block` (offset inside the uncompressed block, below 2**16) to compute the 64-bit virtual offset (`Bio.bgzf.make_virtual_offset()`), or give `virtual_offset` to split it into the two parts (`Bio.bgzf.split_virtual_offset()`). Exactly one of the two forms must be used.
- **API:** `dnakit.sequence.biopython_p18_functions.bgzf_virtual_offset(block_start=None, within_block=None, *, virtual_offset=None)`.
- **Input:** `block_start` and `within_block`, or `virtual_offset`.
- **Output:** `direction` (`make` or `split`), `block_start`, `within_block` and `virtual_offset`.
- **Example code:**

```python
from dnakit.sequence.biopython_p18_functions import bgzf_virtual_offset

made = bgzf_virtual_offset(100, 50)
print(made["direction"], made["virtual_offset"])
split = bgzf_virtual_offset(virtual_offset=made["virtual_offset"])
print(split["direction"], split["block_start"], split["within_block"])
```

- **Example output:**

```text
make 6553650
split 100 50
```
