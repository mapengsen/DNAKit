# Optional molecular-lab backends: DNA design, records and drawing

This group of functions in `dnakit.molecular_lab` is a set of adapters over optional packages: [DnaChisel](https://github.com/Edinburgh-Genome-Foundry/DnaChisel), [DnaCauldron](https://github.com/Edinburgh-Genome-Foundry/DnaCauldron), [pydna](https://github.com/BjornFJohansson/pydna), [pySBOL3](https://github.com/SynBioDex/pySBOL3), [snapgene_reader](https://github.com/Edinburgh-Genome-Foundry/SnapGeneReader), [DNA Features Viewer](https://github.com/Edinburgh-Genome-Foundry/DnaFeaturesViewer), [pyGenomeViz](https://github.com/moshi4/pyGenomeViz) and [primer3-py](https://github.com/libnano/primer3-py). They score sequence constraints, optimize records, design primers and overlap tails, find recombination sites, build and read SBOL3 documents, read SnapGene files, draw feature and genome maps, and compute primer Tm, hairpin and dimer energies. None of them predicts a wet-lab outcome, queries a remote database, or replaces the native Golden Gate and Gibson simulators.

primer3-py is GPL-2.0. DNAKit imports it only inside the `calculate_primer3_*()` functions and copies none of its code.

The dependencies are optional and are imported only when the matching function is called. If a package is missing, the function raises an error that includes the install command. The examples import from `dnakit.molecular_lab`.

## Example data

The only declared input is a one-record FASTA file with a `GGTCTC` site at positions 3 to 9. The other examples pass their inputs directly, or write GenBank and SBOL3 files to a temporary directory.

**Example file:** `part.fa`

```text
>part1
ATGGGTCTCATGCCCAAAGCTTTAA
```

## 1) `dnakit.molecular_lab.evaluate_dnachisel_specs()`

- **Purpose:** Score DnaChisel constraints and objectives without changing the sequence.
- **Method:** Build the requested DnaChisel specifications and evaluate them with `constraints_evaluations()` and `objectives_evaluations()`. Each row reports whether the specification passes, its score, its 0-based locations and DnaChisel's message.
- **API:** `dnakit.molecular_lab.evaluate_dnachisel_specs(sequence, *, constraints=(), objectives=(), circular=False)`.
- **Input:** A DNA string. Each constraint or objective is a mapping with a `name` from the allow-list (`AvoidPattern`, `AvoidChanges`, `AvoidHairpins`, `AvoidStopCodons`, `CodonOptimize`, `EnforceChanges`, `EnforceGCContent`, `EnforceMeltingTemperature`, `EnforceTranslation`, `MaximizeCAI`) and keyword arguments.
- **Output:** `all_constraints_pass`, `failed_constraints`, and the `constraints` and `objectives` rows.
- **Example code:**

```python
from dnakit.molecular_lab import evaluate_dnachisel_specs

result = evaluate_dnachisel_specs(
    "ATGGGTCTCATGCCCAAA",
    constraints=[
        {"name": "AvoidPattern", "pattern": "GGTCTC"},
        {"name": "EnforceGCContent", "mini": 0.3, "maxi": 0.7, "window": 18},
    ],
)
print(result["all_constraints_pass"], result["failed_constraints"])
for row in result["constraints"]:
    print(row["name"], row["passes"], row["locations"])
```

- **Example output:**

```text
False ['AvoidPattern[0-18](pattern:GGTCTC)']
AvoidPattern False [[3, 9]]
EnforceGCContent True []
```

## 2) `dnakit.molecular_lab.optimize_dnachisel_record()`

- **Purpose:** Optimize a local FASTA, GenBank or SnapGene record and write the result as an annotated GenBank file.
- **Method:** Load the record with DnaChisel, resolve the constraints, optimize the objectives when they are given, and write the result with `to_record()`. The output path must be new.
- **API:** `dnakit.molecular_lab.optimize_dnachisel_record(source, *, output_path, constraints=(), objectives=(), circular=False)`.
- **Input:** A local record file, a new `.gb`, `.gbk` or `.genbank` output path, and the same specification mappings as above.
- **Output:** The sequence length, the number of edits, whether all constraints pass, and the output path.
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.molecular_lab import optimize_dnachisel_record

output = Path(tempfile.mkdtemp()) / "part_optimized.gb"
result = optimize_dnachisel_record(
    "part.fa",
    output_path=output,
    constraints=[{"name": "AvoidPattern", "pattern": "GGTCTC"}],
)
print(result["sequence_length"], result["edits"], result["all_constraints_pass"])
```

- **Example output:**

```text
25 1 True
```

## 3) `dnakit.molecular_lab.reverse_translate_protein()`

- **Purpose:** Back-translate a protein into DNA.
- **Method:** Use DnaChisel's `reverse_translate()` with the first codon in table order for each residue. The result is deterministic. A trailing `*` becomes a stop codon. No codon usage is applied.
- **API:** `dnakit.molecular_lab.reverse_translate_protein(protein, *, table="Standard", max_residues=10000)`.
- **Input:** A protein string of the 20 standard residues, optionally ending in `*`, and a codon table name from DnaChisel's list.
- **Output:** `sequence`, `sequence_length`, `protein_length` and `table`.
- **Example code:**

```python
from dnakit.molecular_lab import reverse_translate_protein

result = reverse_translate_protein("MAK*")
print(result["sequence"], result["sequence_length"])
```

- **Example output:**

```text
ATGGCTAAATAA 12
```

## 4) `dnakit.molecular_lab.read_dnacauldron_parts()`

- **Purpose:** Load local part files through DnaCauldron and list the overhangs found in their annotations.
- **Method:** Call DnaCauldron's `load_record()` for each file, then `list_overhangs_from_record_annotations()`. Part IDs must be unique. The total sequence length is bounded.
- **API:** `dnakit.molecular_lab.read_dnacauldron_parts(sources, *, topology="default_to_linear", max_parts=50)`.
- **Input:** A list of GenBank, FASTA or SnapGene file paths.
- **Output:** For each part: `id`, `file`, `length`, `topology`, `overhangs` and `sequence`.
- **Example code:**

```python
from dnakit.molecular_lab import read_dnacauldron_parts

result = read_dnacauldron_parts(["part.fa"])
part = result["parts"][0]
print(part["id"], part["length"], part["topology"])
```

- **Example output:**

```text
part1 25 linear
```

## 5) `dnakit.molecular_lab.write_dnacauldron_constructs()`

- **Purpose:** Simulate a DnaCauldron Gibson or Type IIS assembly and write each construct as an annotated GenBank file.
- **Method:** Build a DnaCauldron repository from the named parts, run `GibsonAssembly` or `Type2sRestrictionAssembly`, and write each construct with `write_record()`. The output directory is created by the call and must not exist. `simulate_dnacauldron_assembly()` in `dnakit.integrations` returns only sequences; this function keeps the feature annotations in the files.
- **API:** `dnakit.molecular_lab.write_dnacauldron_constructs(parts, *, output_dir, method="gibson", enzyme="auto", max_constructs=20)`.
- **Input:** A mapping from part name to DNA, a new output directory, `method` (`"gibson"` or `"type2s"`), and the Type IIS enzyme name when `method="type2s"`.
- **Output:** The `constructs` list with each file name, length and topology, plus any DnaCauldron warnings and errors.
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.molecular_lab import write_dnacauldron_constructs

a = "ATGGCTAAAGGCGAATTGTTTACCGGCGTGGTCCCGATTC"
s1 = "TGGAACTGGATGGCGATGTG"
b = "AACGGCCATAAATTTAGCGTGAGCGGCGAAGGCGAAGGTG"
s2 = "CCGGGTTCAAGTCTCCAACA"
c = "TTAGGCCAACGTTCGATTAGCATCGATCGA"
parts = {"A": a + s1, "B": s1 + b + s2, "C": s2 + c + a[:20]}
result = write_dnacauldron_constructs(
    parts, output_dir=Path(tempfile.mkdtemp()) / "constructs"
)
construct = result["constructs"][0]
print(len(result["constructs"]), construct["length"], construct["topology"])
```

- **Example output:**

```text
1 150 circular
```

## 6) `dnakit.molecular_lab.design_pcr_primers()`

- **Purpose:** Design a forward and a reverse primer that span a linear DNA template.
- **Method:** Call pydna's `primer_design()` on the template. The primers start at the template ends, so the product is the whole template and `product_length` equals the template length.
- **API:** `dnakit.molecular_lab.design_pcr_primers(template, *, target_tm=55.0, min_primer_length=13)`.
- **Input:** A DNA string of at least `2 * min_primer_length` bases. `target_tm` is a temperature between 40 and 75 °C.
- **Output:** `forward_primer`, `reverse_primer`, their lengths and `product_length`.
- **Example code:**

```python
from dnakit.molecular_lab import design_pcr_primers

template = "ATGGCTAAAGGCGAATTGTTTACCGGCGTGGTCCCGATTCTGGTGGAACTGGATGGCGATGTGAACGG"
result = design_pcr_primers(template)
print(result["forward_primer"], result["reverse_primer"], result["product_length"])
```

- **Example output:**

```text
ATGGCTAAAGGCGAATTG CCGTTCACATCGCC 68
```

## 7) `dnakit.molecular_lab.design_assembly_primers()`

- **Purpose:** Add overlap tails to the primers of several PCR fragments so that adjacent fragments can be joined by Gibson assembly or fusion PCR.
- **Method:** Design a primer pair for each template with pydna's `primer_design()`, then call pydna's `assembly_fragments()`. The reverse complement of fragment i's reverse primer ends with the first `overlap` bases of fragment i+1's forward primer.
- **API:** `dnakit.molecular_lab.design_assembly_primers(templates, *, overlap=35, max_link=40, circular=False)`.
- **Input:** Two to 20 DNA templates, an overlap of 10 to 100 bases, and `max_link` of at least `overlap`.
- **Output:** For each fragment: `forward_primer`, `reverse_primer` and `product_length`.
- **Example code:**

```python
from dnakit.molecular_lab import design_assembly_primers

fragment_1 = (
    "ATGGCTAAAGGCGAATTGTTTACCGGCGTGGTCCCGATTCTGGTGGAACTGGATGGCGATGTGAACGGCCATAAATTTAGC"
)
fragment_2 = (
    "GTGAGCGGCGAAGGCGAAGGTGATGCACCGTACCCGGAACTGACCCTGAAATTCATTTGCACCACCGGCAAACTGCCGGTG"
)
result = design_assembly_primers([fragment_1, fragment_2], overlap=20)
for item in result["fragments"]:
    print(item["index"], item["forward_primer"], item["reverse_primer"], item["product_length"])
```

- **Example output:**

```text
0 ATGGCTAAAGGCGAATTG CGCCGCTCACGCTAAATTTATGGCCGTTCA 91
1 TAAATTTAGCGTGAGCGGCGAAG CACCGGCAGTTTGC 91
```

## 8) `dnakit.molecular_lab.find_recombination_sites()`

- **Purpose:** Find loxP or Gateway att sites in a DNA sequence.
- **Method:** Call pydna's `find_loxP_sites()` for `system="loxp"`, or `find_gateway_sites()` for `system="gateway"`. The result is sorted by start position.
- **API:** `dnakit.molecular_lab.find_recombination_sites(sequence, *, system="loxp", greedy=False)`.
- **Input:** A DNA string. `system` is `"loxp"` or `"gateway"`. `greedy` is passed to the Gateway finder.
- **Output:** `site_count` and `sites`. Each site has `name`, `start`, `end` (0-based, end-exclusive) and `strand`.
- **Example code:**

```python
from dnakit.molecular_lab import find_recombination_sites

loxp = find_recombination_sites("GGG" + "ATAACTTCGTATAGCATACATTATACGAAGTTAT" + "CCC")
print(loxp["site_count"], loxp["sites"])
attb = find_recombination_sites("AAACAAGTTTGTACAAAAAAGCAGGCTCCCC", system="gateway")
print(attb["site_count"], attb["sites"])
```

- **Example output:**

```text
1 [{'name': 'loxP', 'start': 3, 'end': 37, 'strand': 1}]
1 [{'name': 'attB1', 'start': 3, 'end': 24, 'strand': 1}]
```

## 9) `dnakit.molecular_lab.build_sbol3_document()`

- **Purpose:** Build an SBOL3 document in which each named DNA sequence is a Component with one Sequence.
- **Method:** Create a pySBOL3 `Sequence` with IUPAC DNA encoding and a DNA `Component` for each display ID. Identities are `namespace + display_id` and `namespace + display_id + "_sequence"`. The document is validated before it is serialized.
- **API:** `dnakit.molecular_lab.build_sbol3_document(components, *, namespace, file_format="turtle", output_path=None)`.
- **Input:** A mapping from valid SBOL display IDs to DNA strings, an `http(s)` namespace ending in `/` or `#`, and `file_format` set to `turtle`, `ntriples`, `rdfxml` or `jsonld`.
- **Output:** `component_count`, `sequence_bases`, `validation` (error and warning counts with details), and `document` (the text, when `output_path` is not given) or `output_path`.
- **Example code:**

```python
from dnakit.molecular_lab import build_sbol3_document

result = build_sbol3_document(
    {"promoter_1": "TTGACAGCTAGCTCAGTCCTAGGTATAATGCTAGC", "cds_1": "ATGGCTAAAGGCGAATAA"},
    namespace="https://dnakit.example.org/parts/",
)
print(result["component_count"], result["sequence_bases"], result["validation"]["error_count"])
```

- **Example output:**

```text
2 53 0
```

## 10) `dnakit.molecular_lab.read_sbol3_sequences()`

- **Purpose:** Read the DNA sequences stored in a local SBOL3 document.
- **Method:** Load the file with pySBOL3's `Document.read()` and return every top-level `Sequence` object with its identity, display ID, encoding and elements.
- **API:** `dnakit.molecular_lab.read_sbol3_sequences(source, *, max_objects=50000, max_sequences=1000)`.
- **Input:** A local SBOL3 file in a format pySBOL3 can infer from its extension.
- **Output:** `object_count`, `sequence_count` and `sequences`.
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.molecular_lab import build_sbol3_document, read_sbol3_sequences

path = Path(tempfile.mkdtemp()) / "parts.ttl"
build_sbol3_document(
    {"promoter_1": "TTGACAGCTAGC", "cds_1": "ATGGCTAAATAA"},
    namespace="https://dnakit.example.org/parts/",
    output_path=path,
)
for row in read_sbol3_sequences(path)["sequences"]:
    print(row["display_id"], row["length"], row["elements"])
```

- **Example output:**

```text
cds_1_sequence 12 ATGGCTAAATAA
promoter_1_sequence 12 TTGACAGCTAGC
```

## 11) `dnakit.molecular_lab.validate_sbol3_document()`

- **Purpose:** Validate a local SBOL3 document and list its errors and warnings.
- **Method:** Load the file with pySBOL3 and call `Document.validate()`. Each issue reports the object, the rule and the message.
- **API:** `dnakit.molecular_lab.validate_sbol3_document(source, *, max_issues=200)`.
- **Input:** A local SBOL3 file.
- **Output:** `valid`, `error_count`, `warning_count`, and up to `max_issues` entries in `errors` and `warnings`.
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.molecular_lab import build_sbol3_document, validate_sbol3_document

path = Path(tempfile.mkdtemp()) / "parts.ttl"
build_sbol3_document(
    {"promoter_1": "TTGACAGCTAGC"},
    namespace="https://dnakit.example.org/parts/",
    output_path=path,
)
checked = validate_sbol3_document(path)
print(checked["valid"], checked["error_count"], checked["warning_count"])
```

- **Example output:**

```text
True 0 0
```

## 12) `dnakit.molecular_lab.read_snapgene_sequence()`

- **Purpose:** Read the sequence, topology and methylation flags from a SnapGene `.dna` file.
- **Method:** Call snapgene_reader's `snapgene_file_to_dict()` and map its fields to DNAKit names. snapgene_reader only reads files. DNAKit does not write `.dna` files.
- **API:** `dnakit.molecular_lab.read_snapgene_sequence(source, *, max_bases=2000000)`.
- **Input:** A SnapGene `.dna` file written by SnapGene.
- **Output:** `sequence`, `length`, `topology`, `strandedness`, the three methylation flags and the export and import versions.
- **Example code (not run locally):**

```python
from dnakit.molecular_lab import read_snapgene_sequence

result = read_snapgene_sequence("pUC19.dna")
print(result["length"], result["topology"], result["sequence"][:20])
```

- **Example output (not run locally):**

```text
<length> <topology> <first 20 bases>
```

The example is not run locally because no SnapGene-written `.dna` file is available to the checker. The output block is a template, not a recorded result.

## 13) `dnakit.molecular_lab.read_snapgene_features()`

- **Purpose:** List the features of a SnapGene `.dna` file with their coordinates.
- **Method:** Read the `features` block from snapgene_reader and convert each start to a 0-based, end-exclusive interval. The label comes from the `label` qualifier, or from the feature name when that qualifier is missing.
- **API:** `dnakit.molecular_lab.read_snapgene_features(source, *, max_features=10000)`.
- **Input:** A SnapGene `.dna` file.
- **Output:** `feature_count` and `features`, each with `name`, `type`, `start`, `end`, `strand`, `color` and `label`.
- **Example code (not run locally):**

```python
from dnakit.molecular_lab import read_snapgene_features

result = read_snapgene_features("pUC19.dna")
for row in result["features"][:3]:
    print(row["label"], row["start"], row["end"], row["strand"])
```

- **Example output (not run locally):**

```text
<label> <start> <end> <strand>
```

The example is not run locally, for the same reason as the previous function.

## 14) `dnakit.molecular_lab.convert_snapgene_to_genbank()`

- **Purpose:** Convert a SnapGene `.dna` file to a new GenBank file.
- **Method:** Call snapgene_reader's `snapgene_file_to_gbk()`. The GenBank text comes from that converter, and its LOCUS date is fixed by it. On failure, the partial output is removed.
- **API:** `dnakit.molecular_lab.convert_snapgene_to_genbank(source, output_path)`.
- **Input:** A SnapGene `.dna` file and a new `.gb`, `.gbk` or `.genbank` output path.
- **Output:** `output_path` and the number of bytes written.
- **Example code (not run locally):**

```python
from dnakit.molecular_lab import convert_snapgene_to_genbank

result = convert_snapgene_to_genbank("pUC19.dna", "pUC19.gb")
print(result["output_path"], result["bytes"])
```

- **Example output (not run locally):**

```text
<output path> <bytes>
```

The example is not run locally, for the same reason as the previous function.

## 15) `dnakit.molecular_lab.build_feature_map()`

- **Purpose:** Return feature-map data for a linear sequence: coordinates, colors and the lane each feature is drawn on.
- **Method:** Create DNA Features Viewer `GraphicFeature` objects and call `compute_features_levels()`. Overlapping features are moved to higher lanes. Lanes start at 0.0.
- **API:** `dnakit.molecular_lab.build_feature_map(sequence_length, features, *, max_features=10000)`.
- **Input:** The sequence length and a list of mappings with `start` and `end` (0-based, end-exclusive), and optional `strand` (1, -1 or 0), `label` and `color` (`#RRGGBB`).
- **Output:** `lane_count` and `features`, each with its coordinates, strand, color, label and `lane`.
- **Example code:**

```python
from dnakit.molecular_lab import build_feature_map

result = build_feature_map(
    100,
    [
        {"start": 0, "end": 40, "label": "promoter"},
        {"start": 30, "end": 60, "strand": -1, "label": "tag"},
        {"start": 70, "end": 90, "label": "terminator"},
    ],
)
for row in result["features"]:
    print(row["label"], row["start"], row["end"], row["lane"])
```

- **Example output:**

```text
promoter 0 40 0.0
tag 30 60 1.0
terminator 70 90 0.0
```

## 16) `dnakit.molecular_lab.render_feature_map_svg()`

- **Purpose:** Draw a linear feature map and write it as a new SVG file.
- **Method:** Draw a DNA Features Viewer `GraphicRecord` on a Matplotlib figure with an Agg canvas, so no display is needed. The SVG has no date metadata, and its element IDs use a fixed hash salt, so the same input gives the same bytes.
- **API:** `dnakit.molecular_lab.render_feature_map_svg(sequence_length, features, *, output_path, figure_width=8.0, max_features=10000)`.
- **Input:** The same features as `build_feature_map()`, plus a new `.svg` output path. `figure_width` is between 2 and 20 inches.
- **Output:** `output_path`, `format`, `bytes`, `feature_count` and `sequence_length`.
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.molecular_lab import render_feature_map_svg

path = Path(tempfile.mkdtemp()) / "map.svg"
result = render_feature_map_svg(
    100,
    [
        {"start": 0, "end": 40, "label": "promoter"},
        {"start": 30, "end": 60, "strand": -1, "label": "tag"},
    ],
    output_path=path,
)
print(result["format"], result["feature_count"], path.read_bytes()[:5])
```

- **Example output:**

```text
svg 2 b'<?xml'
```

## 17) `dnakit.molecular_lab.build_genbank_feature_map()`

- **Purpose:** Return the full feature-map data for a local GenBank record, with the lane of every feature.
- **Method:** Read the record with DNA Features Viewer's `BiopythonTranslator` and compute lanes with `compute_features_levels()`. `summarize_dna_features_viewer()` in `dnakit.integrations` returns only a preview of the first features.
- **API:** `dnakit.molecular_lab.build_genbank_feature_map(source, *, max_features=10000)`.
- **Input:** A local GenBank file.
- **Output:** `sequence_length`, `feature_count`, `lane_count` and `features`, each with `label`, coordinates, `strand`, `color` and `lane`.
- **Example code:**

```python
import tempfile
from pathlib import Path

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import FeatureLocation, SeqFeature
from Bio.SeqRecord import SeqRecord

from dnakit.molecular_lab import build_genbank_feature_map

record = SeqRecord(Seq("A" * 100), id="map1", name="map1", description="example")
record.annotations["molecule_type"] = "DNA"
record.features = [
    SeqFeature(FeatureLocation(10, 40, strand=1), type="CDS", qualifiers={"gene": ["geneX"]}),
    SeqFeature(FeatureLocation(30, 60, strand=-1), type="gene", qualifiers={"gene": ["geneY"]}),
]
path = Path(tempfile.mkdtemp()) / "map.gb"
SeqIO.write(record, str(path), "genbank")
result = build_genbank_feature_map(path)
print(result["sequence_length"], result["lane_count"])
print([(row["label"], row["start"], row["end"]) for row in result["features"]])
```

- **Example output:**

```text
100 2
[('geneX', 10, 40), ('geneY', 30, 60)]
```

## 18) `dnakit.molecular_lab.read_genbank_genome_tracks()`

- **Purpose:** Read the sequence length and selected features of each GenBank record as a genome track.
- **Method:** Get the sizes with pyGenomeViz's `Genbank` parser, then select features of one key from the parsed records. The filter is DNAKit's own loop, because pyGenomeViz 0.4.4's `get_seqid2features()` reads `SeqFeature.strand`, which Biopython 1.88 no longer provides. Multi-part locations use their overall span.
- **API:** `dnakit.molecular_lab.read_genbank_genome_tracks(source, *, feature_type="CDS", max_records=100, max_features=50000)`.
- **Input:** A local GenBank file and one feature key, such as `"CDS"` or `"gene"`.
- **Output:** `track_count`, `feature_count` and `tracks`. Each track has `name`, `length` and `features`, and each feature has `type`, `start`, `end`, `strand` and `label`.
- **Example code:**

```python
import tempfile
from pathlib import Path

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import FeatureLocation, SeqFeature
from Bio.SeqRecord import SeqRecord

from dnakit.molecular_lab import read_genbank_genome_tracks

record = SeqRecord(Seq("A" * 100), id="map1", name="map1", description="example")
record.annotations["molecule_type"] = "DNA"
record.features = [
    SeqFeature(FeatureLocation(10, 40, strand=1), type="CDS", qualifiers={"gene": ["geneX"]}),
    SeqFeature(FeatureLocation(30, 60, strand=-1), type="gene", qualifiers={"gene": ["geneY"]}),
]
path = Path(tempfile.mkdtemp()) / "map.gb"
SeqIO.write(record, str(path), "genbank")
result = read_genbank_genome_tracks(path, feature_type="CDS")
track = result["tracks"][0]
print(result["track_count"], track["length"], track["features"])
```

- **Example output:**

```text
1 100 [{'type': 'CDS', 'start': 10, 'end': 40, 'strand': 1, 'label': 'geneX'}]
```

## 19) `dnakit.molecular_lab.render_genome_tracks_svg()`

- **Purpose:** Draw feature tracks, with optional links between adjacent tracks, and write them as a new SVG file.
- **Method:** Add one pyGenomeViz feature track per entry, add its features, then add the links with `add_link()`. Links are allowed only between neighbouring tracks in the order given, because pyGenomeViz draws links only there. pyGenomeViz writes a creation timestamp and random group IDs into the SVG, so the bytes differ between runs.
- **API:** `dnakit.molecular_lab.render_genome_tracks_svg(tracks, *, output_path, links=(), fig_width=15.0, max_tracks=20)`.
- **Input:** A list of tracks, each with `name`, `length` and `features`. The optional `links` are mappings with `track_a`, `start_a`, `end_a`, `track_b`, `start_b` and `end_b`.
- **Output:** `output_path`, `format`, `bytes`, `track_count`, `feature_count` and `link_count`.
- **Example code:**

```python
import tempfile
from pathlib import Path

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import FeatureLocation, SeqFeature
from Bio.SeqRecord import SeqRecord

from dnakit.molecular_lab import read_genbank_genome_tracks, render_genome_tracks_svg

record = SeqRecord(Seq("A" * 100), id="map1", name="map1", description="example")
record.annotations["molecule_type"] = "DNA"
record.features = [
    SeqFeature(FeatureLocation(10, 40, strand=1), type="CDS", qualifiers={"gene": ["geneX"]}),
]
path = Path(tempfile.mkdtemp()) / "map.gb"
SeqIO.write(record, str(path), "genbank")
track = read_genbank_genome_tracks(path)["tracks"][0]
other = {"name": "other", "length": 120, "features": [{"start": 5, "end": 60, "strand": -1}]}
result = render_genome_tracks_svg(
    [track, other],
    output_path=Path(tempfile.mkdtemp()) / "tracks.svg",
    links=[
        {
            "track_a": track["name"],
            "start_a": 10,
            "end_a": 40,
            "track_b": "other",
            "start_b": 5,
            "end_b": 60,
        }
    ],
)
print(result["format"], result["track_count"], result["feature_count"], result["link_count"])
```

- **Example output:**

```text
svg 2 2 1
```

## 20) `dnakit.molecular_lab.calculate_primer3_tm()`

- **Purpose:** Compute the melting temperature of one primer with primer3-py.
- **Method:** Call primer3-py's `calc_tm()` with the SantaLucia nearest-neighbour model. Monovalent and divalent salt, dNTP and primer concentrations are passed through. primer3-py is GPL-2.0; DNAKit imports it only when this function runs.
- **API:** `dnakit.molecular_lab.calculate_primer3_tm(sequence, *, mv_conc=50.0, dv_conc=1.5, dntp_conc=0.6, dna_conc=50.0)`.
- **Input:** A primer of 1 to 200 bases using A, C, G and T only. Concentrations are in mM (monovalent, divalent and dNTP) and nM (primer).
- **Output:** `tm_c` in degrees Celsius, `sequence_length` and the `conditions` used.
- **Example code:**

```python
from dnakit.molecular_lab import calculate_primer3_tm

result = calculate_primer3_tm("ACGTTGCAAGCTTGCA")
print(result["sequence_length"], round(result["tm_c"], 2))
```

- **Example output:**

```text
16 54.92
```

## 21) `dnakit.molecular_lab.calculate_primer3_hairpin()`

- **Purpose:** Report the hairpin stability of one primer with primer3-py.
- **Method:** Call primer3-py's `calc_hairpin()`. The thermodynamic fields are None when primer3 finds no hairpin. Energies are reported in cal/mol.
- **API:** `dnakit.molecular_lab.calculate_primer3_hairpin(sequence, *, mv_conc=50.0, dv_conc=1.5, dntp_conc=0.6, dna_conc=50.0, temp_c=37.0, max_loop=30)`.
- **Input:** A primer as in `calculate_primer3_tm()`. `temp_c` is between 0 and 100 °C, and `max_loop` is between 3 and 100.
- **Output:** `structure_found`, `tm_c`, `dg_cal_per_mol`, `dh_cal_per_mol` and `ds_cal_per_mol_k`.
- **Example code:**

```python
from dnakit.molecular_lab import calculate_primer3_hairpin

result = calculate_primer3_hairpin("GGGGAAAACCCC")
print(result["structure_found"], round(result["tm_c"], 2), round(result["dg_cal_per_mol"], 1))
```

- **Example output:**

```text
True 55.6 -1312.3
```

## 22) `dnakit.molecular_lab.calculate_primer3_dimer()`

- **Purpose:** Report the self-dimer of one primer, or the hetero-dimer with a second primer, using primer3-py.
- **Method:** Call primer3-py's `calc_homodimer()` when no partner is given, and `calc_heterodimer()` otherwise. The thermodynamic fields are None when primer3 finds no dimer.
- **API:** `dnakit.molecular_lab.calculate_primer3_dimer(sequence, *, partner=None, mv_conc=50.0, dv_conc=1.5, dntp_conc=0.6, dna_conc=50.0, temp_c=37.0, max_loop=30)`.
- **Input:** A primer, an optional `partner` primer, and the same conditions as `calculate_primer3_hairpin()`.
- **Output:** `dimer_type` (`"homodimer"` or `"heterodimer"`), `structure_found`, `tm_c`, `dg_cal_per_mol`, `dh_cal_per_mol`, `ds_cal_per_mol_k` and `partner_length`.
- **Example code:**

```python
from dnakit.molecular_lab import calculate_primer3_dimer

result = calculate_primer3_dimer("ATGCGTACGT", partner="ACGTACGCAT")
print(result["dimer_type"], result["structure_found"], round(result["tm_c"], 2))
```

- **Example output:**

```text
heterodimer True 31.61
```

## Installing the optional dependencies

```bash
python -m pip install "dnakit[bio-dnachisel]"            # dnachisel: functions 1 to 3
python -m pip install "dnakit[bio-dnacauldron]"          # dnacauldron: functions 4 and 5
python -m pip install "dnakit[bio-assembly]"             # pydna: functions 6 to 8
python -m pip install "dnakit[bio-pysbol3]"              # sbol3: functions 9 to 11
python -m pip install "dnakit[bio-snapgene]"             # snapgene_reader: functions 12 to 14
python -m pip install "dnakit[bio-dna-features-viewer]" "dnakit[bio-plots]"  # functions 15 to 17
python -m pip install "dnakit[bio-pygenomeviz]"          # pyGenomeViz: functions 18 and 19
python -m pip install "dnakit[bio-primer3]"              # primer3-py: functions 20 to 22
```

Project pages: [DnaChisel](https://github.com/Edinburgh-Genome-Foundry/DnaChisel), [DnaCauldron](https://github.com/Edinburgh-Genome-Foundry/DnaCauldron), [pydna](https://github.com/BjornFJohansson/pydna), [pySBOL3](https://github.com/SynBioDex/pySBOL3), [snapgene_reader](https://github.com/Edinburgh-Genome-Foundry/SnapGeneReader), [DNA Features Viewer](https://github.com/Edinburgh-Genome-Foundry/DnaFeaturesViewer), [pyGenomeViz](https://github.com/moshi4/pyGenomeViz), [primer3-py](https://github.com/libnano/primer3-py).

## How these functions are verified

The test file is `others/tests/unit/test_molecular_lab_backends.py`. It has 41 items: 40 run and one is skipped unless `DNAKIT_SNAPGENE_TEST_FILE` names a SnapGene-written `.dna` file. The checks are:

- Hand calculations: the DnaChisel pattern location `[3, 9)`, the first-codon back-translation `ATGGCTAAATAA`, the loxP site at `[3, 37)`, the lanes `[0, 1, 0]` for three features, the 370 bp circular Gibson product, the 303 bp GenBank track, and the overlap of 20 bases between adjacent primers from `design_assembly_primers()`.
- Direct upstream calls: DnaChisel constraint evaluation, pydna primer design and Gateway site search, DNA Features Viewer levels, pyGenomeViz's GenBank size parser, and the primer3-py Tm, hairpin and dimer functions.
- Parsing and validation through the real upstream: SBOL3 round trips and an invalid untyped Component, and rejection of non-SnapGene bytes by snapgene_reader. Bad sequences, bad namespaces, non-adjacent pyGenomeViz links and existing output paths are refused.
- Determinism: two renders of the SVG feature map give identical bytes.
- A missing optional package raises `BackendUnavailableError` with the extra name, for each of the eight packages.

## Current scope

These functions wrap the DNA-level operations that the existing adapters did not cover: constraint scoring and record output for DnaChisel; file input and annotated GenBank output for DnaCauldron; primer design, overlap tails and site search for pydna; SBOL3 construction, sequence reading and validation; SnapGene sequence and feature reading and conversion; feature-map data and SVG output for DNA Features Viewer; GenBank tracks and SVG output for pyGenomeViz; and primer3 Tm, hairpin and dimer values.

Not wrapped, and why:

- DnaChisel optimization from a sequence string (`optimize_with_dnachisel()`), DnaCauldron simulation from strings (`simulate_dnacauldron_assembly()`), SBOL3 inventory (`summarize_sbol3_document()`) and feature-map previews (`summarize_dna_features_viewer()`) already exist in `dnakit.integrations`.
- pydna's Gibson product (`assemble_duplexes()`) already exists in `dnakit.bioinformatics.cloning`.
- primer3-py's `designPrimers()`: primer design already exists through `dnakit.molbio.Primer3CLIDesignAdapter`, which calls the `primer3_core` program. This batch wraps only the thermodynamic functions that are missing.
- SnapGene writing: snapgene_reader is read-only, so no `.dna` file is written.
- DnaChisel's BLAST helper (it needs network access), pydna's fusion-PCR helper and gel simulation (`simulate_gel()` already exists in `dnakit.bioinformatics.cloning`), DnaCauldron's PDF reports and mix planning, and pyGenomeViz's HTML export and GC or skew profiles. The last group needs an interactive or extra output, or is left for a later batch.
- The SnapGene read examples are not run locally, because the repository holds no SnapGene-written `.dna` file. The reader's error path is tested with a non-SnapGene file.
