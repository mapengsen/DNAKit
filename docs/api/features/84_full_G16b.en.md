# 84. Third-party package interfaces G16b: SBOL3, SnapGene, DNA Features Viewer, pyGenomeViz, primer3-py

This page wraps the DNA- or gene-level public functions of five upstream packages that run offline, all under `dnakit.molecular_lab`. There are 12 wrapped functions. Each returns a JSON-serializable dictionary with `status`, `operation`, the result fields and `provenance`.

- Optional dependencies: `bio-pysbol3`, `bio-snapgene`, `bio-dna-features-viewer`, `bio-pygenomeviz` and `bio-primer3`. Upstream packages are imported only inside the functions; a missing package raises `BackendUnavailableError`.
- License: primer3-py is GPL-2.0. It is imported only optionally, and dnakit does not copy its code.
- Network: `fetch_genbank_by_accession()` calls NCBI. Its tests replace the client with a fake, so the default test run makes no network request.
- SnapGene: no `.dna` sample file exists locally, so the positive-path test is skipped; the reason is stated here.
- Not wrapped: plotting and layout functions, command-line entry points, example-data downloads, workflows that need external programs (MMseqs2, MUMmer, progressiveMauve), deprecated aliases, and functions already covered by existing dnakit interfaces (`backends.py` or `dnakit.ops`).

## copy_sbol3_design() {#84-sbol3-001}

- **Purpose:** Copy the top-level SBOL3 objects of a local file into a new namespace, optionally writing the copies to a new file.
- **Method:** Clones the top-level objects with `sbol3.copy`, keeping display IDs and paths. `sbol3.copy(into_namespace=...)` can keep a reference to the old namespace when a referencing object is processed before the object it references; this function repeats the reference update over the whole group, so the result does not depend on input order.
- **API:** `dnakit.molecular_lab.copy_sbol3_design(source, *, into_namespace, output=None, verbose=True)`
- **Input:** Path to an SBOL file (`.ttl`, `.nt`, `.rdf`, `.xml` or `.jsonld`). `into_namespace` must be an absolute URI. `output` is a new file path; its suffix selects the serialization.
- **Output:** `copied_count` and the `copied` list sorted by identity (type, source identity, new identity, display ID), plus `output`, `output_format` and `provenance`.
- **Example file:** `design.ttl`

```turtle
@prefix sbol: <http://sbols.org/v3#> .

<https://example.org/a/part1> a sbol:Component ;
    sbol:displayId "part1" ;
    sbol:hasNamespace <https://example.org/a> ;
    sbol:hasSequence <https://example.org/a/seq1> ;
    sbol:type <https://identifiers.org/SBO:0000251> .

<https://example.org/a/seq1> a sbol:Sequence ;
    sbol:displayId "seq1" ;
    sbol:elements "ATGCGTACGATCG" ;
    sbol:encoding <https://identifiers.org/edam:format_1207> ;
    sbol:hasNamespace <https://example.org/a> .
```

- **Example code:**

```python
import pathlib

from dnakit.molecular_lab import copy_sbol3_design

pathlib.Path("copied.nt").unlink(missing_ok=True)
result = copy_sbol3_design(
    "design.ttl", into_namespace="https://example.org/b", output="copied.nt"
)
print(result["copied_count"], result["output_format"])
for row in result["copied"]:
    print(row["type"], row["identity"])
```

- **Example output:**

```text
2 NTRIPLES
Component https://example.org/b/part1
Sequence https://example.org/b/seq1
```

## to_sbol3_display_id() {#84-sbol3-002}

- **Purpose:** Convert free text into a valid SBOL3 displayId, for naming DNA parts.
- **Method:** Calls `sbol3.string_to_display_id`, which replaces disallowed characters and prefixes an underscore when the text starts with a digit; the result is checked with `sbol3.is_valid_display_id`. Blank text raises an error.
- **API:** `dnakit.molecular_lab.to_sbol3_display_id(name, *, verbose=True)`
- **Input:** A non-empty string of at most 256 characters.
- **Output:** `display_id`, `valid`, `rule` and `provenance`.

- **Example code:**

```python
from dnakit.molecular_lab import to_sbol3_display_id

print(to_sbol3_display_id("pTac promoter 1")["display_id"])
result = to_sbol3_display_id("1abc")
print(result["display_id"], result["valid"])
```

- **Example output:**

```text
pTac_promoter_1
_1abc True
```

## check_sbol3_display_id() {#84-sbol3-003}

- **Purpose:** Check whether a string is a valid SBOL3 displayId.
- **Method:** Calls `sbol3.is_valid_display_id`. Rule: letters, digits and underscores only, not starting with a digit.
- **API:** `dnakit.molecular_lab.check_sbol3_display_id(display_id, *, verbose=True)`
- **Input:** A string; an empty string is reported as invalid.
- **Output:** `display_id`, `valid`, `rule` and `provenance`.

- **Example code:**

```python
from dnakit.molecular_lab import check_sbol3_display_id

for text in ["pTac_1", "1bad", "a-b"]:
    print(text, check_sbol3_display_id(text)["valid"])
```

- **Example output:**

```text
pTac_1 True
1bad False
a-b False
```

## extract_sbol3_display_id() {#84-sbol3-004}

- **Purpose:** Read the displayId from an SBOL3 URL identity, which is the last segment of the URL path.
- **Method:** Calls `sbol3.extract_display_id`. A non-URL identity returns `None`; a last segment that is not a valid displayId raises `ConfigurationError`.
- **API:** `dnakit.molecular_lab.extract_sbol3_display_id(identity, *, verbose=True)`
- **Input:** The identity string of an SBOL object.
- **Output:** `identity`, `display_id`, `found` and `provenance`.

- **Example code:**

```python
from dnakit.molecular_lab import extract_sbol3_display_id

print(extract_sbol3_display_id("https://example.org/lab/pTac")["display_id"])
print(extract_sbol3_display_id("plain text")["found"])
```

- **Example output:**

```text
pTac
False
```

## snapgene_seqrecord_summary() {#84-snapgene-001}

- **Purpose:** Read a SnapGene `.dna` file with snapgene_reader and return a summary of the Biopython SeqRecord view (sequence and features).
- **Method:** Calls `snapgene_file_to_seqrecord`. Sequence length and feature count are bounded by `max_bases` and `max_features`; features are sorted by position. The file is only read. A missing or unparseable file raises `ConfigurationError`.
- **API:** `dnakit.molecular_lab.snapgene_seqrecord_summary(source, *, max_bases=2000000, max_features=10000, verbose=True)`
- **Input:** Path to a local `.dna` file.
- **Output:** `record_id`, `name`, `length`, `sequence`, `feature_count` and `features` sorted by position (type, start, end, strand, label).
- **Note:** No SnapGene `.dna` sample exists in the local environment, so the positive-path test is skipped; the missing-file and malformed-file tests are kept. The example below shows only the error path.

- **Example code:**

```python
import pathlib

from dnakit.exceptions import ConfigurationError
from dnakit.molecular_lab import snapgene_seqrecord_summary

pathlib.Path("bad.dna").write_bytes(b"not a SnapGene file")
for name in ["missing.dna", "bad.dna"]:
    try:
        snapgene_seqrecord_summary(name)
    except ConfigurationError as exc:
        print(name, exc.code, exc.message.split(":")[0])
```

- **Example output:**

```text
missing.dna CONFIGURATION_ERROR Input file does not exist
bad.dna CONFIGURATION_ERROR snapgene_reader could not parse bad.dna
```

## annotate_sequence_feature() {#84-dfv-001}

- **Purpose:** Add one annotated feature to a DNA sequence and return that feature.
- **Method:** Builds a Biopython SeqRecord and calls `dna_features_viewer.annotate_biopython_record`. Positions are 0-based, with start inclusive and end exclusive.
- **API:** `dnakit.molecular_lab.annotate_sequence_feature(sequence, *, start=0, end=None, strand=1, feature_type="misc_feature", qualifiers=None, verbose=True)`
- **Input:** A DNA string (IUPAC ambiguity codes allowed). `end` defaults to the sequence end; `strand` is 1, -1 or 0; `qualifiers` is a string-to-string mapping.
- **Output:** `sequence_length` and `feature` (type, start, end, strand, qualifiers), plus `provenance`.

- **Example code:**

```python
from dnakit.molecular_lab import annotate_sequence_feature

result = annotate_sequence_feature(
    "ATGCCCGGGTAA", start=3, end=9, feature_type="CDS", qualifiers={"label": "core"}
)
feature = result["feature"]
print(result["sequence_length"], feature["type"], feature["start"], feature["end"])
print(feature["strand"], feature["qualifiers"]["label"])
```

- **Example output:**

```text
12 CDS 3 9
1 core
```

## graphical_translation() {#84-dfv-002}

- **Purpose:** Translate a DNA subsegment in the left-to-right order used by the feature map.
- **Method:** Calls `dna_features_viewer.biotools.extract_graphical_translation`. On the reverse strand the reverse complement is translated and the result is reversed, so the amino acids follow the forward sequence.
- **API:** `dnakit.molecular_lab.graphical_translation(sequence, *, start=0, end=None, strand=1, long_form=False, verbose=True)`
- **Input:** A DNA string; `strand` is 1 or -1; `long_form=True` returns a list of three-letter residue names.
- **Output:** `translation`, `residue_count`, `start`, `end`, `strand`, `long_form` and `provenance`.

- **Example code:**

```python
from dnakit.molecular_lab import graphical_translation

print(graphical_translation("ATGGCCTAA")["translation"])
print(graphical_translation("TTAGCCAT", end=8, strand=-1)["translation"])
print(graphical_translation("TTAGCCATG", strand=-1, long_form=True)["translation"])
```

- **Example output:**

```text
MA*
AM
['*', 'Gly', 'His']
```

## summarize_genbank_records() {#84-pgv-001}

- **Purpose:** Read one or more GenBank files with pyGenomeViz's parser and summarize genome length, range, GC percentage and records.
- **Method:** Calls `pygenomeviz.scripts.gbk_files2gbk_objects`. Each input may be written as `path:start-end` (a range) or `path::-1` (reverse complement). Records are sorted by sequence ID.
- **API:** `dnakit.molecular_lab.summarize_genbank_records(sources, *, verbose=True)`
- **Input:** A list of 1 to 100 GenBank file paths.
- **Output:** `genome_count` and `genomes` (name, full length, range, range size, reverse flag, GC percentage, records), plus `provenance`.
- **Example file:** `genome.gb`

```text
LOCUS       TEST1                     49 bp    DNA              UNK 01-JAN-1980
DEFINITION  Example genome fixture.
ACCESSION   TEST1
VERSION     TEST1.1
KEYWORDS    .
SOURCE      .
  ORGANISM  .
            .
FEATURES             Location/Qualifiers
     CDS             1..12
ORIGIN
        1 atgcgtacga tcgtttaagg ccatgcaagt cgatcgatcg ggatcctag
//
```

- **Example code:**

```python
from dnakit.molecular_lab import summarize_genbank_records

result = summarize_genbank_records(["genome.gb", "genome.gb:10-30"])
for genome in result["genomes"]:
    print(genome["name"], genome["full_length"], genome["range"], genome["gc_percent"])
```

- **Example output:**

```text
genome 49 [0, 49] 51.020408
genome 49 [10, 30] 45.0
```

## fetch_genbank_by_accession() {#84-pgv-002}

- **Purpose:** Fetch one nucleotide GenBank record from NCBI by accession. This is a network function.
- **Method:** Calls `pygenomeviz.utils.fetch_genbank_by_accid` (NCBI Entrez efetch) and parses the returned text as one Biopython record. The accession and email are checked before any request; a failed request raises `BackendExecutionError`. Tests replace the Entrez client with a fake, so no live request is made.
- **API:** `dnakit.molecular_lab.fetch_genbank_by_accession(accession, *, output=None, email=None, verbose=True)`
- **Input:** An accession such as `NC_000913.3` or `JX128258.1`; `output` is a new file path; `email` is an optional contact address.
- **Output:** `accession`, `record_id`, `name`, `length`, `output`, `bytes` and `provenance`.

- **Example code:**

```python
from dnakit.exceptions import ConfigurationError
from dnakit.molecular_lab import fetch_genbank_by_accession

try:
    fetch_genbank_by_accession("bad accession")
except ConfigurationError as exc:
    print(exc.code, exc.message)
```

- **Example output:**

```text
CONFIGURATION_ERROR accession must look like NC_000913.3, JX128258.1 or CP000001.
```

## calculate_primer3_end_stability() {#84-p3-001}

- **Purpose:** Compute the 3' end stability of a primer annealed to a template.
- **Method:** Calls `primer3.calc_end_stability` (SantaLucia nearest-neighbour model). Returns Tm, dG, dH and dS; the thermodynamic values are `None` when no structure is found. primer3-py is GPL-2.0 and is only imported optionally.
- **API:** `dnakit.molecular_lab.calculate_primer3_end_stability(primer, template, *, mv_conc=50.0, dv_conc=1.5, dntp_conc=0.6, dna_conc=50.0, temp_c=37.0, max_loop=30, verbose=True)`
- **Input:** Primer and template of 1 to 200 bases (A, C, G, T, N). `mv_conc`, `dv_conc` and `dntp_conc` are in mM; `dna_conc` is in nM.
- **Output:** `structure_found`, `tm_c`, `dg_cal_per_mol`, `dh_cal_per_mol`, `ds_cal_per_mol_k`, `conditions` and `provenance`.

- **Example code:**

```python
from dnakit.molecular_lab import calculate_primer3_end_stability

result = calculate_primer3_end_stability("GCTAGCTCAGTCCTAGG", "CCTAGGACTGAGCTAGC")
print(result["structure_found"], round(result["tm_c"], 2), round(result["dg_cal_per_mol"], 1))
```

- **Example output:**

```text
True 52.24 -17299.9
```

## design_primer3_pairs() {#84-p3-002}

- **Purpose:** Design PCR primer pairs for a DNA template with primer3.
- **Method:** Calls `primer3.design_primers`. An optional `target=(start, length)` (0-based) marks a region the amplicon must cover. Pairs are returned in primer3 penalty order, lowest first. The 5' start of the reverse primer is derived from its 3' position and length.
- **API:** `dnakit.molecular_lab.design_primer3_pairs(template, *, target=None, min_product=100, max_product=300, num_return=5, min_size=18, opt_size=20, max_size=25, min_tm=57.0, opt_tm=60.0, max_tm=63.0, min_gc=40.0, max_gc=60.0, verbose=True)`
- **Input:** A template of 50 to 10000 bp (A, C, G, T, N).
- **Output:** `pair_count` and `pairs`. Each pair holds the left and right oligo sequences, start, length, Tm, GC percentage, product size and penalty.
- **Note:** The example template contains a self-complementary segment, so the left and right oligos of pair 2 are identical. primer3 still returns the pair in penalty order; check primer-dimer and self-complementarity before use.

- **Example code:**

```python
from dnakit.molecular_lab import design_primer3_pairs

template = "ATGCGTACGATCGATCGTAGCTAGCTAGCATCGATCGTAGCTACGATCGTAGCTAGCTAGCTAGCTGATCGATCGATCGTAGCTAGCTAGCATGCATGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCATCGATGCATGCATGCTAGCTAGCTAGCATCGATCGATGCTAGCTAGCTAGCTAGCATGCTAGCTAGCTAGCTACGATCGATCGATCGTAGCTACGATCGA"
result = design_primer3_pairs(template, min_product=80, max_product=200, num_return=2)
print(result["pair_count"])
for pair in result["pairs"]:
    print(pair["rank"], pair["product_size"], pair["left"]["sequence"], pair["right"]["sequence"])
```

- **Example output:**

```text
2
1 189 AGCTAGCATCGATCGTAGCT CGATCGATCGTAGCTAGCTAGC
2 142 CGATCGATCGTAGCTAGCTAGC CGATCGATCGTAGCTAGCTAGC
```

## sanitize_primer3_sequence() {#84-p3-003}

- **Purpose:** Replace IUPAC ambiguity codes with N so that a sequence meets primer3's input rules.
- **Method:** Calls `primer3.p3helpers.sanitize_sequence`. A, C, G, T and N are kept as they are, including case. Characters that are not IUPAC codes, such as X, Z or gaps, raise an error.
- **API:** `dnakit.molecular_lab.sanitize_primer3_sequence(sequence, *, verbose=True)`
- **Input:** A non-empty string.
- **Output:** `input_length`, `sanitized`, `changed_count` and `provenance`.

- **Example code:**

```python
from dnakit.molecular_lab import sanitize_primer3_sequence

result = sanitize_primer3_sequence("ACGTRYN")
print(result["sanitized"], result["changed_count"])
```

- **Example output:**

```text
ACGTNNN 2
```
