# Sequence concatenation and assembly

Splice multiple DNA sequences or perform de-overlapping splicing based on precise overlapping regions.

## 1) Ordinary sequence splicing

- **Purpose:** Connect multiple DNA fragments in the input order, insert linker or structured gap at the junction, and return the coordinates of each fragment in the new sequence, suitable for constructing combined sequences.
- **Method:** Join fragments in input order, inserting the configured linker or Gap at each junction.
- **API:** `dnakit.ops.concat(sequences[required], linker[optional], gap[optional])`.
- **Input:** At least two `DNASequence` or normalized strings; `linker` and `gap` cannot be used at the same time.
- **Output:** Returns the concatenated `DNASequence`.
- **Example code:**

```python
from dnakit import DNASequence
from dnakit.ops import concat

joined = concat(
    [DNASequence("AA"), DNASequence("CC")],
    linker="T",
)
print(joined.symbols)
```

- **Example output:**

```text
AATCC
```

## 2) Remove overlapping splicing

- **Purpose:** Complete the connection based on the overlapping relationship between the ends of two pieces of DNA, so that the matched overlapping region is retained only once, and the joint position and final sequence are reported.
- **Method:** Choose the longest exact left-suffix/right-prefix overlap and retain one copy; compare IUPAC symbols literally.
- **API:** `dnakit.ops.concat_overlap(sequences[required], min_overlap[optional], max_overlap[optional])`.
- **Input:** Must be two linear `DNASequence`s with no explicit gap or a normalized string. `min_overlap=1` by default; `max_overlap` limits the overlap search.
- **Output:** Returns the processed `DNASequence`; read its symbols from `symbols`.
- **Example code:**

```python
from dnakit import DNASequence
from dnakit.ops import concat_overlap

joined = concat_overlap(
    [DNASequence("AAACCC"), DNASequence("CCCGG")],
)
print(joined.symbols)
```

- **Example output:**

```text
AAACCCGG
```

<!-- BEGIN GENERATED API: cloning -->

<span id="cloning"></span>

## 3) Design overlaps for Gibson assembly {#dnakit.databases.dna_gibson_design}

- **Purpose:** Design overlaps for Gibson assembly.
- **Method:** Take overlaps of the requested length from neighboring fragment ends and append them to assembly parts.
- **API：** `dnakit.databases.dna_gibson_design(*, operation=None, fragments, overlap_length=None, verbose=True)`
- **Input:** Provide `fragments`.
- **Output:** Returns a mapping with `assembly_fragments`, `n_fragments`, `overlap_length` and related fields.
- **Example code:**

```python
import dnakit

example_fragments = ['ATGGAGGAGCCGCAGTCAGATCCTAGCGTTGAATTCGGATCCAAGCTTACGTACGT',
 'CGTACGTAAGCTTGGATCCGAATTCAAGCGATAGGATCTGACTGCGGCTCCTCCAT',
 'TCCATGGAGGAGCCGCAGTCAGATCCTAGCGTTGAATTCGGATCCAAGCTTACGTACGT']

result = dnakit.databases.dna_gibson_design(
    operation='gibson_design',
    fragments=example_fragments,
    overlap_length=20,
    verbose=False,
)
print(result)
```

- **Example output:** The response has the following shape; records depend on the input and source.

Response-shape excerpt; `<…>` denotes type placeholders, not a captured query result.

```json
{
  "assembly_fragments": [
    {
      "name": "<string>",
      "original_sequence": "<string>",
      "with_overlaps": "<string>",
      "left_overlap": "<string>",
      "right_overlap": "<string>"
    }
  ],
  "n_fragments": "<integer>",
  "overlap_length": "<integer>"
}
```

[Full parameters and usage notes](../../faq/api/04_ops_006_sequence_concatenation.md#dnakit.databases.dna_gibson_design-parameters)

## 4) Determine Golden Gate reaction products {#dnakit.databases.dna_golden_gate_assemble}

- **Purpose:** Determine Golden Gate reaction products.
- **Method:** Simulate Type IIS digestion, discard fragments retaining recognition sites and ligate matching overhangs.
- **API：** `dnakit.databases.dna_golden_gate_assemble(*, operation=None, fragments, enzyme=None, circular=None, labels=None, verbose=True)`
- **Input:** Provide `fragments`.
- **Output:** Returns a mapping with `enzyme`, `product_length`, `num_fragments_assembled`, `assembly_order` and related fields.
- **Example code:**

```python
import dnakit

example_fragments = ['GGTCTCAAAACATGGAGGAGCCGCAGTCAGATCCTAGCGTTGAATTCGGATCCCTTTTGAGACC',
 'GGTCTCAAAAGGGATCCAAGCTTACGTACGTATGGAGGAGCCGCAGTCAGATATTTTGAGACC',
 'GGTCTCAAAATCCTAGCGTTGAATTCGGATCCAAGCTTACGTACGTATGGAGCGTTTGAGACC',
 'GGTCTCAAACGATCGATCGATCGATCGATCGGTTTTGAGACC']

result = dnakit.databases.dna_golden_gate_assemble(
    fragments=example_fragments,
    enzyme='BsaI',
    circular=True,
    verbose=False,
)
print(result)
```

- **Example output:** The response has the following shape; records depend on the input and source.

Response-shape excerpt; `<…>` denotes type placeholders, not a captured query result.

```json
{
  "enzyme": "<string>",
  "product_length": "<integer>",
  "num_fragments_assembled": "<integer>",
  "assembly_order": [
    "<value>"
  ],
  "unused_fragments": [
    "<value>"
  ],
  "product_sequence": "<string>"
}
```

[Full parameters and usage notes](../../faq/api/04_ops_006_sequence_concatenation.md#dnakit.databases.dna_golden_gate_assemble-parameters)

## 5) Design Golden Gate restriction sites and overhangs {#dnakit.databases.dna_golden_gate_design}

- **Purpose:** Design Golden Gate restriction sites and overhangs.
- **Method:** Assign unique nonpalindromic four-base overhangs and add Type IIS recognition sites.
- **API：** `dnakit.databases.dna_golden_gate_design(*, operation=None, parts, enzyme=None, verbose=True)`
- **Input:** Provide `parts`.
- **Output:** Returns a mapping with `parts_with_overhangs`, `overhangs`, `enzyme`, `protocol_note` and related fields.
- **Example code:**

```python
import dnakit

example_parts = ['ATGGAGGAGCCGCAGTCAGATCCTAGCGTTGAATTCGGATCC',
 'GGATCCAAGCTTACGTACGTATGGAGGAGCCGCAGTCAGAT',
 'CCTAGCGTTGAATTCGGATCCAAGCTTACGTACGTATGGAG']

result = dnakit.databases.dna_golden_gate_design(
    operation='golden_gate_design',
    parts=example_parts,
    enzyme='BsaI',
    verbose=False,
)
print(result)
```

- **Example output:** The response has the following shape; records depend on the input and source.

Response-shape excerpt; `<…>` denotes type placeholders, not a captured query result.

```json
{
  "parts_with_overhangs": [
    {
      "name": "<string>",
      "sequence": "<string>",
      "left_overhang": "<string>",
      "right_overhang": "<string>",
      "full_sequence": "<string>"
    }
  ],
  "overhangs": [
    "<string>"
  ],
  "enzyme": "<string>",
  "protocol_note": "<string>"
}
```

[Full parameters and usage notes](../../faq/api/04_ops_006_sequence_concatenation.md#dnakit.databases.dna_golden_gate_design-parameters)

<!-- END GENERATED API: cloning -->

## 6) Assemble reads with SPAdes {#bio-18}

- **Purpose:** Assemble short reads or hybrid short/long reads into contigs.
- **Method:** Run SPAdes; `second` supplies mate reads and `nanopore` supplies hybrid long reads.
- **API:** `dnakit.bioinformatics.ngs.assemble_reads(first, output_directory, second=None, nanopore=None, memory_gb=8)`
- **Input:** FASTQ input, SPAdes on PATH, a new output directory and sufficient memory for `memory_gb`.
- **Output:** Returns `ToolResult` with `backend`, artifact paths in `artifacts`, and execution `log`.
- **Agent/MCP API:** `dnakit.ngs.assemble_reads(read1, output_directory, read2=None, nanopore=None, memory_gb=8, ...)`; returns JSON-serializable assembly artifact paths and diagnostics.
- **Citation and license:** cite SPAdes as Bankevich et al. (2012), DOI `10.1089/cmb.2012.0021` (GPL-2.0). The caller installs the backend.
- **Example code:**

```python
from dnakit.bioinformatics.ngs import assemble_reads

result = assemble_reads("R1.fq.gz", "new_assembly", second="R2.fq.gz", memory_gb=8)
print(sorted(result.artifacts))
```

- **Example output:**

Returns `ToolResult`; `artifacts` lists SPAdes outputs such as contigs and scaffolds, with keys determined by actual outputs. This external-program example was not executed.

[Dependencies and implementation notes](../../faq/bioinformatics.md#reads).

## 7) Ligate duplex DNA fragments {#bio-52-ligate}

- **Purpose:** Join duplex DNA fragments with compatible ends.
- **Method:** Check end polarity and complementarity, then join both strands.
- **API:** `dnakit.bioinformatics.cloning.Duplex.from_sequence(sequence, circular=False)`; `Duplex.ligate(other)`
- **Input:** Two linear `Duplex` objects with compatible ends; this example cuts one molecule and ligates its fragments.
- **Output:** Returns the ligated `Duplex` with its Watson strand and end information.
- **Example code:**

```python
from dnakit.bioinformatics.cloning import Duplex

molecule = Duplex.from_sequence("AACCGGTT")
left, right = molecule.cut(3, 5)
print(left.ligate(right).watson)
```

- **Example output:**

```text
AACCGGTT
```

[Dependencies and implementation notes](../../faq/bioinformatics.md#optimization).

## 8) Record cloning steps and intermediate products {#bio-53}

- **Purpose:** Organize cloning operations into traceable named steps.
- **Method:** Register products with `add`; use `apply` to run functions and store parents and products.
- **API:** `dnakit.bioinformatics.cloning.CloningPlan()`; `CloningPlan.add(name, product)`; `CloningPlan.apply(name, function, parents=(), **parameters)`
- **Input:** Python callables, parent-step names and operation parameters.
- **Output:** Returns a `CloningPlan`; the example shows `add`, `apply`, `steps`.
- **Example code:**

```python
from dnakit.bioinformatics.cloning import CloningPlan, Duplex

plan = CloningPlan()
plan.add("template", Duplex.from_sequence("AACCGGTT"))
cut = plan.apply("cut", lambda molecule: molecule.cut(3, 5), parents=["template"])
print([step.name for step in plan.steps])
print(len(cut.product))
```

- **Example output:**

```text
['template', 'cut']
2
```

[Dependencies and implementation notes](../../faq/bioinformatics.md#optimization).

## 9) Assemble overlapping duplex fragments {#bio-54}

- **Purpose:** Find linear or circular products assembled from overlapping fragments.
- **Method:** Use pydna to search assembly paths satisfying the minimum overlap length.
- **API:** `dnakit.bioinformatics.cloning.assemble_duplexes(fragments, circular=False, minimum_overlap=25, max_assemblies=50)`
- **Input:** A list of `Duplex` fragments; requires `bio-assembly`; set overlap length according to the experimental design.
- **Output:** Returns a tuple of assembled `Duplex` products.
- **Example code:**

```python
from dnakit.bioinformatics.cloning import Duplex, assemble_duplexes

fragments = [Duplex.from_sequence("TTATACGATGCTAGCTTACGGTCA"),
             Duplex.from_sequence("AGCTTACGGTCACCTATGAGTCCA")]
products = assemble_duplexes(fragments, minimum_overlap=12)
longest = max(products, key=lambda product: len(product.watson))
print(longest.watson)
```

- **Example output:**

```text
TTATACGATGCTAGCTTACGGTCACCTATGAGTCCA
```

[Dependencies and implementation notes](../../faq/bioinformatics.md#optimization).

## 10) Golden Gate part design {#13-golden-gate-part-design}

- **Purpose:** Assign unique non-palindromic junction overhangs and add BsaI/BbsI sites.
- **Method:** Assign unique non-palindromic overhangs and add BsaI/BbsI recognition sites.
- **API:** `dnakit.molbio.design_golden_gate(parts[required], enzyme[optional])`
- **Input:** At least two canonical DNA values; design supports `BsaI` and `BbsI`.
- **Output:** Returns `ProviderResult` with query content in `data`, plus `provider` and provenance.
- **Example code:**

```python
from dnakit import DNASequence
from dnakit.molbio import design_golden_gate

result = design_golden_gate(
    (DNASequence("ATGGCTGAA"), DNASequence("GCCAAATAA")),
    enzyme="BsaI",
)
print(result.data["parts_with_overhangs"])
```

- **Example output:** Left/right overhangs, full synthetic part sequences, and a protocol note.

[Optional backend dependencies and return conventions](../../faq/tutorial_integrations.md#optional-bioinformatics)

## 11) Golden Gate reaction assembly {#14-golden-gate-reaction-assembly}

- **Purpose:** Simulate Type IIS digestion and ligation by matching overhangs.
- **Method:** Simulate Type IIS digestion and join fragments by compatible overhangs.
- **API:** `dnakit.molbio.assemble_golden_gate(fragments[required], enzyme[optional], circular[optional], labels[optional])`
- **Input:** At least two canonical DNA values carrying the selected Type IIS sites; BsaI, BbsI, Esp3I/BsmBI, and SapI are supported.
- **Output:** Returns `ProviderResult` with query content in `data`, plus `provider` and provenance.
- **Example code:**

```python
from dnakit import DNASequence
from dnakit.molbio import assemble_golden_gate

result = assemble_golden_gate(
    (
        DNASequence("GGTCTCAAAACATGGAGGAGCCGCAGTCAGATCCTAGCGTTGAATTCGGATCCCTTTTGAGACC"),
        DNASequence("GGTCTCAAAAGGGATCCAAGCTTACGTACGTATGGAGGAGCCGCAGTCAGATATTTTGAGACC"),
        DNASequence("GGTCTCAAAATCCTAGCGTTGAATTCGGATCCAAGCTTACGTACGTATGGAGCGTTTGAGACC"),
        DNASequence("GGTCTCAAACGATCGATCGATCGATCGATCGGTTTTGAGACC"),
    ),
    labels=("vector", "insert-1", "insert-2", "insert-3"),
)
print(result.data["product_sequence"])
```

- **Example output:** Product sequence, length, assembly order, junction overhangs, and unused fragments.

[Optional backend dependencies and return conventions](../../faq/tutorial_integrations.md#optional-bioinformatics)
