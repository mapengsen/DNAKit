# G18n review wrappers: pydna, DnaCauldron and DnaChisel

These review wrappers cover the functions that the ledgers still listed as `other` for pydna, DnaCauldron and DnaChisel. Each wrapper sits in one of three modules under `src/dnakit/molecular_lab/`: `pydna_review_functions.py`, `dnacauldron_review_functions.py` and `dnachisel_review_functions.py`. The lead wires them into the package later, so the examples import from the module files; after wiring the same names are exported from `dnakit.molecular_lab`.

Each function imports its package only when it runs. A missing package raises `BackendUnavailableError` with the install hint. Every function returns a JSON-safe dictionary with `status`, `operation`, the result fields and `provenance`. Coordinates are 0-based; intervals are half-open. All examples run offline. Functions that stay unwrapped, and why, are listed at the end.

## 1) `dnakit.molecular_lab.dna_molecules_equal()`

- **Purpose:** Test whether DNA strings describe the same double-stranded molecule.
- **Method:** Calls pydna's `eq` on the strings. A sequence equals its reverse complement; with `circular=True` every rotation also counts as equal. Upstream: `pydna.utils.eq`.
- **API:** `dnakit.molecular_lab.dna_molecules_equal(sequences, *, circular=False, verbose=True)`
- **Input:** At least two DNA strings (A, C, G, T; whitespace removed, case ignored) and a boolean `circular`.
- **Output:** `equal`, `sequence_count`, `lengths`, `circular` and `provenance`.
- **Example code:**

```python
from dnakit.molecular_lab import dna_molecules_equal
result = dna_molecules_equal(["ATGC", "GCAT"])
print(result["equal"], dna_molecules_equal(["ATGC", "TGCA"])["equal"])
print(dna_molecules_equal(["ATGCAAGTCC", "AAGTCCATGC"], circular=True)["equal"])
```

- **Example output:**

```text
True False
True
```

## 2) `dnakit.molecular_lab.rare_codon_positions()`

- **Purpose:** Locate the codons that are rare in an organism's codon usage table.
- **Method:** Calls pydna's `rarecodons` with the bundled `sce` (yeast) or `eco` (E. coli) rare-codon table. Codons are read in frame from position 0; trailing bases that do not form a codon are ignored.
- **API:** `dnakit.molecular_lab.rare_codon_positions(sequence, *, organism="sce", verbose=True)`
- **Input:** A DNA string and `organism`, either `sce` or `eco`.
- **Output:** `codon_count`, `rare_codon_count`, `positions` (each with `start`, `end` and `codon`) and `provenance`.
- **Example code:**

```python
from dnakit.molecular_lab import rare_codon_positions
result = rare_codon_positions("ATGCGATAGAGAAAAAAAAGCTCGA")
print(result["codon_count"], result["rare_codon_count"])
print(result["positions"])
```

- **Example output:**

```text
8 1
[{'start': 3, 'end': 6, 'codon': 'CGA'}]
```

## 3) `dnakit.molecular_lab.smallest_circular_rotation()`

- **Purpose:** Return the canonical form of a circular DNA sequence: its lexicographically smallest rotation.
- **Method:** Calls pydna's `smallest_rotation`, which uses a Lyndon-word algorithm (pydivsufsort is used when installed).
- **API:** `dnakit.molecular_lab.smallest_circular_rotation(sequence, *, verbose=True)`
- **Input:** A DNA string.
- **Output:** `sequence` (the rotated sequence), `length` and `provenance`.
- **Example code:**

```python
from dnakit.molecular_lab import smallest_circular_rotation
print(smallest_circular_rotation("CAGT")["sequence"])
print(smallest_circular_rotation("GATTACA")["sequence"])
```

- **Example output:**

```text
AGTC
ACAGATT
```

## 4) `dnakit.molecular_lab.generate_random_orf()`

- **Purpose:** Generate a random open reading frame: ATG, random sense codons and a stop codon.
- **Method:** Calls pydna's `randomORF`. The seed is set for this call only; the global random state is saved before and restored after the call. With `max_codons`, the length is drawn from a triangular distribution between the two values.
- **API:** `dnakit.molecular_lab.generate_random_orf(codons, *, max_codons=None, seed=0, verbose=True)`
- **Input:** `codons`, an integer from 2 to 10000 counting the start and stop codons; optional `max_codons`, at least `codons`.
- **Output:** `sequence`, `length` (3 times `codons`), `codon_count`, `seed`, `max_codons` and `provenance`.
- **Example code:**

```python
from dnakit.molecular_lab import generate_random_orf
result = generate_random_orf(10, seed=1)
print(result["sequence"], result["codon_count"])
```

- **Example output:**

```text
ATGACGGACGCAGTGTCTCTGTCGATATAG 10
```

## 5) `dnakit.molecular_lab.annotate_primer_binding_sites()`

- **Purpose:** Mark where two primers bind the two ends of a product.
- **Method:** Checks that the forward primer starts the product and that the reverse primer, written 5' to 3', is the reverse complement of the product's end. Then calls pydna's `annotate_primer_binding_sites`, which adds a `primer_bind` feature at each end.
- **API:** `dnakit.molecular_lab.annotate_primer_binding_sites(product, forward_primer, reverse_primer, *, verbose=True)`
- **Input:** DNA strings for the product and both primers. A mismatch raises `ConfigurationError`.
- **Output:** `product_length`, `feature_count`, `features` (each with `type`, `start`, `end`, `strand` and `label`) and `provenance`.
- **Example code:**

```python
from dnakit.molecular_lab import annotate_primer_binding_sites
result = annotate_primer_binding_sites(
    "ATGCGTACGTTAGCATCGGATCCAAAGGGCCC", "ATGCGTACGTTAGCAT", "GGGCCCTTT"
)
for feature in result["features"]:
    print(feature["label"], feature["start"], feature["end"], feature["strand"])
```

- **Example output:**

```text
forward_primer 0 16 1
reverse_primer 23 32 -1
```

## 6) `dnakit.molecular_lab.suggest_pcr_program()`

- **Purpose:** Simulate PCR with a primer pair and return pydna's suggested cycling programme.
- **Method:** Runs pydna's `amplify.pcr`, where each primer must anneal over at least `min_annealing` bases. Then applies the Tm and annealing rule of the chosen polymerase: `taq` (`program`), `dbd` (`dbd_program`) or `pfu_sso7d` (`pfu_sso7d_program`).
- **API:** `dnakit.molecular_lab.suggest_pcr_program(template, forward_primer, reverse_primer, *, polymerase="taq", min_annealing=14, verbose=True)`
- **Input:** Template DNA, two primers, `polymerase` and `min_annealing` from 8 to 60.
- **Output:** `product_length`, `program` (the text block), `program_lines`, `polymerase` and `provenance`.
- **Example code:**

```python
from dnakit.molecular_lab import suggest_pcr_program
A = "ATGGCTAAAGGTGAAGAACTTTTCACTGGAGTTGTCCCAATTCTTGTTGAATTAGATGGTGATGTTAATGGGCACAAATTTTCTGTCAGTGGAGAGGGTGAAGG"
result = suggest_pcr_program(A, A[:22], A[-22:].translate(str.maketrans("ACGT", "TGCA"))[::-1])
print(result["product_length"], result["polymerase"])
for line in result["program_lines"]:
    print(line)
```

- **Example output:**

```text
104 taq
|95°C|95°C               |    |tmf:58.1
|____|_____          72°C|72°C|tmr:65.4
|3min|30s  \ 51.8°C _____|____|45s/kb
|    |      \______/ 0:30|5min|GC 41%
|    |       30s         |    |104bp
```

## 7) `dnakit.molecular_lab.ligate_blunt_fragments()`

- **Purpose:** Enumerate the products of ligating blunt, double-stranded DNA fragments.
- **Method:** Calls pydna's `ligation_assembly` with `allow_blunt=True`. Sticky ends from restriction digests are not represented. Products are sorted by length and then sequence; `product_count` gives the total.
- **API:** `dnakit.molecular_lab.ligate_blunt_fragments(fragments, *, allow_partial_overlap=False, circular_only=False, max_products=50, verbose=True)`
- **Input:** 2 to 20 blunt DNA strings; `max_products` from 1 to 200.
- **Output:** `fragment_count`, `product_count`, `products` (each with `sequence`, `length` and `circular`) and `provenance`.
- **Example code:**

```python
from dnakit.molecular_lab import ligate_blunt_fragments
first = "ATGCGTACGTTAGCATCGGATCCAAAGGGCCCTTTAAAGGAGATCGTTTAGCGCCC"
second = "TTTAGCGCCCAAAGCTTGCATGCCTGCAGGTCGACTCTAGAGGATCCCCGGGTACCGAGC"
result = ligate_blunt_fragments([first, second])
print(result["product_count"], [row["length"] for row in result["products"]])
print(all(row["circular"] for row in result["products"]))
```

- **Example output:**

```text
2 [116, 116]
True
```

## 8) `dnakit.molecular_lab.gateway_recombination_products()`

- **Purpose:** Enumerate Gateway BP or LR recombination products of DNA fragments.
- **Method:** Calls pydna's `gateway_assembly` for the chosen reaction. `greedy=True` uses the wider site definitions; otherwise the conservative sites are used. Fragments must carry the att sites in their DNA.
- **API:** `dnakit.molecular_lab.gateway_recombination_products(fragments, *, reaction_type="BP", greedy=False, circular_only=False, multi_site_only=False, max_products=50, verbose=True)`
- **Input:** 2 to 20 DNA strings with att sites; `reaction_type` is `BP` or `LR`.
- **Output:** `reaction_type`, `site_set`, `product_count`, `products` and `provenance`.
- **Example code:**

```python
from dnakit.molecular_lab import gateway_recombination_products
from Bio.Data.IUPACData import ambiguous_dna_values
from pydna.gateway import recombinase_dict
site = recombinase_dict["BP"]["conservative"].recombinases[0]
def expand(text):
    return "".join(ambiguous_dna_values[ch][0] if ch in ambiguous_dna_values else ch for ch in text.upper())
att_b, att_p = expand(site.site1), expand(site.site2)
linear = att_b + "TTTTACGTACGTAAAGGGCCCAAATTTGGGGCCAAGGCCTTAAGGCCAA" + att_b
plasmid = att_p + "AAAACCCGGGTTTAAACCCGGG" + att_p
result = gateway_recombination_products([linear, plasmid], reaction_type="BP")
print(result["reaction_type"], result["site_set"], result["product_count"])
print([row["length"] for row in result["products"]])
```

- **Example output:**

```text
BP conservative 7
[96, 157, 166, 227, 324, 350, 411]
```

## 9) `dnakit.molecular_lab.simulate_genome_edit()`

- **Purpose:** Simulate one genome-editing event on a DNA sequence.
- **Method:** Dispatches to pydna's `homologous_recombination_*`, `cre_lox_*`, `recombinase_*` or `crispr_integration`. Integration needs inserts; excision and inversion take none. CRISPR needs guide spacers that cut the genome. Recombinase needs two sites (upper-case flanks, lower-case homology core) such as `("ATGCCCTAAaaTT", "AAaaTTTTTTTCCCT")`. CRE-lox integration needs circular inserts, so set `inserts_circular=True`. `limit` (minimum homology) applies to homologous recombination and CRISPR only.
- **API:** `dnakit.molecular_lab.simulate_genome_edit(genome, *, mechanism, operation, inserts=(), guides=(), recombinase_sites=None, genome_circular=False, inserts_circular=False, limit=40, max_products=50, verbose=True)`
- **Input:** Genome DNA up to 20000 bases; `mechanism` (`homologous_recombination`, `cre_lox`, `recombinase` or `crispr`); `operation` (`integration`, `excision` or `excision_or_inversion`; CRISPR allows only integration).
- **Output:** `mechanism`, `edit`, `genome_length`, `insert_count`, `guide_count`, `limit`, `product_count`, `products` and `provenance`.
- **Example code:**

```python
from dnakit.molecular_lab import simulate_genome_edit
result = simulate_genome_edit(
    "CCCCCCATAACTTCGTATAGCATACATTATACGAAGTTATGGGGGATAACTTCGTATAGCATACATTATACGAAGTTATAAAAA",
    mechanism="cre_lox",
    operation="excision",
)
print(result["edit"], result["product_count"], [row["length"] for row in result["products"]])
```

- **Example output:**

```text
excision 2 [39, 45]
```

- **Example code:**

```python
from dnakit.molecular_lab import simulate_genome_edit
result = simulate_genome_edit(
    "aaccggttcaatgcaaacagtaatgatggatgacattcaaagcac",
    mechanism="crispr",
    operation="integration",
    inserts=["aaccggttAAAAAAAAAttcaaagcac"],
    guides=["ttcaatgcaaacagtaatga"],
    limit=8,
)
print(result["mechanism"], result["limit"], [row["sequence"] for row in result["products"]])
```

- **Example output:**

```text
crispr 8 ['AACCGGTTAAAAAAAAATTCAAAGCAC']
```

## 10) `dnakit.molecular_lab.assemble_recombinase_fragments()`

- **Purpose:** Enumerate recombinase assembly products of DNA fragments.
- **Method:** Calls pydna's `recombinase_assembly` with a two-site recombinase built from `recombinase_sites`, in the same format as `simulate_genome_edit`. Each fragment carries a site in its DNA; matching sites join the fragments into products.
- **API:** `dnakit.molecular_lab.assemble_recombinase_fragments(fragments, *, recombinase_sites, circular_only=False, max_products=50, verbose=True)`
- **Input:** 2 to 20 DNA strings, each carrying a site, and the site pair.
- **Output:** `fragment_count`, `product_count`, `products` and `provenance`.
- **Example code:**

```python
from dnakit.molecular_lab import assemble_recombinase_fragments
sites = ("ATGCCCTAAaaTT", "AAaaTTTTTTTCCCT")
result = assemble_recombinase_fragments(
    ["CCCCATGCCCTAAaaTTGGGGGGGG", "AAaaTTTTTTTCCCTTTTTTTTT"], recombinase_sites=sites
)
print(result["product_count"], [row["length"] for row in result["products"]])
```

- **Example output:**

```text
2 [14, 34]
```

## 11) `dnakit.molecular_lab.extract_protospacers()`

- **Purpose:** Extract Cas9 protospacers that sit directly before the guide scaffold.
- **Method:** Calls pydna's `protospacer` with the `cas9` class: a 20-base spacer followed by the scaffold `GTTTTAGAGCTAGAAATAGCAAGTTAAAATAAGG`. Both strands are searched; `circular=True` also searches across the origin.
- **API:** `dnakit.molecular_lab.extract_protospacers(guide_construct, *, cas="cas9", circular=False, verbose=True)`
- **Input:** A guide construct DNA string; `cas` must be `cas9`, the only nuclease pydna defines.
- **Output:** `protospacer_count`, `protospacers` (sorted), `protospacer_length`, `scaffold`, `circular` and `provenance`.
- **Example code:**

```python
from dnakit.molecular_lab import extract_protospacers
from pydna.crispr import cas9
construct = "ATTA" + "GACGTACGATCGTAGCTAGG" + cas9.scaffold + "TTTT"
result = extract_protospacers(construct)
print(result["protospacer_count"], result["protospacers"])
```

- **Example output:**

```text
1 ['GACGTACGATCGTAGCTAGG']
```

## 12) `dnakit.molecular_lab.parse_sequence_text()`

- **Purpose:** Parse EMBL, GenBank or FASTA text into DNA records.
- **Method:** Calls pydna's `parsers.embl_gb_fasta` on the text in memory, so no file is needed. Records keep their input order.
- **API:** `dnakit.molecular_lab.parse_sequence_text(text, *, verbose=True)`
- **Input:** Text of at most 1 MB that holds at least one record.
- **Output:** `record_count`, `records` (each with `name`, `length`, `topology` and `sequence`) and `provenance`.
- **Example code:**

```python
from dnakit.molecular_lab import parse_sequence_text
result = parse_sequence_text(">one\nATGCCC\n>two\nGGGAAA\n")
for record in result["records"]:
    print(record["name"], record["length"], record["topology"], record["sequence"])
```

- **Example output:**

```text
one 6 linear ATGCCC
two 6 linear GGGAAA
```

## 13) `dnakit.molecular_lab.design_user_cloning_primers()`

- **Purpose:** Design USER-cloning primer tails that give neighbouring fragments overlapping ends.
- **Method:** Amplifies each template with its primers through pydna's `amplify.pcr`, then calls pydna's `user_assembly_design` with `max_overlap`, `min_overlap` and `max_tail`. Uracil in the tails is kept as `U`, as pydna writes it.
- **API:** `dnakit.molecular_lab.design_user_cloning_primers(fragments, *, max_overlap=15, min_overlap=4, max_tail=50, verbose=True)`
- **Input:** 2 to 20 mappings, each with `template`, `forward_primer` and `reverse_primer` DNA strings; the primers must amplify the template.
- **Output:** `fragment_count`, `fragments` (each with `forward_primer`, `reverse_primer` and `product_length`) and `provenance`.
- **Example code:**

```python
from dnakit.molecular_lab import design_user_cloning_primers
A = "ATGGCTAAAGGTGAAGAACTTTTCACTGGAGTTGTCCCAATTCTTGTTGAATTAGATGGTGATGTTAATGGGCACAAATTTTCTGTCAGTGGAGAGGGTGAAGG"
B = "GGCGAAGGGATGCTTACCGTTGGTACACACCCCAATCTTTTACTGGAGCGTAAAGGTTACGCATGGCAATTAAATCGCTGGATTGTCGATCCACGCCATCTGG"
def reverse_complement(text):
    return text.translate(str.maketrans("ACGT", "TGCA"))[::-1]
fragments = [
    {"template": A, "forward_primer": A[:22], "reverse_primer": reverse_complement(A[-22:])},
    {"template": B, "forward_primer": B[:22], "reverse_primer": reverse_complement(B[-22:])},
]
result = design_user_cloning_primers(fragments)
for row in result["fragments"]:
    print(row["forward_primer"], row["reverse_primer"], row["product_length"])
```

- **Example output:**

```text
ATGGCTAAAGGTGAAGAAC ATCCCTTCGCCCCUTCACCCTCTCCA 104
AGGGGCGAAGGGAUGCT CCAGATGGCGTGGAT 103
```

## 14) `dnakit.molecular_lab.simulate_type2s_assembly()`

- **Purpose:** Simulate a Golden Gate (Type-2S) assembly of linear parts into circular products.
- **Method:** Builds DnaCauldron's `generate_type2s_restriction_mix` and asks it for the circular assemblies. `enzyme="auto"` lets DnaCauldron's `autoselect_enzyme` choose among BbsI, AarI, BsaI, BsmBI and SapI.
- **API:** `dnakit.molecular_lab.simulate_type2s_assembly(parts, *, enzyme="auto", max_products=50, verbose=True)`
- **Input:** A mapping from part name to a linear DNA string whose ends carry the overhangs.
- **Output:** `enzyme`, `part_names`, `circular_product_count`, `products` (each with `sequence`, `length` and `circular`) and `provenance`.
- **Example code:**

```python
from dnakit.molecular_lab import simulate_type2s_assembly
INS1 = "ATGGCTAAAGGTGAAGAACTTTTCACTGGAGTTGTCCCAATTCTTGTTGAATTAGATGGTGATGTTAATGGG"
INS2 = "ATGGCTAGCGTGATCAAGCCGGAAATGCGCATCCACCTGCGTAAAAGTGTTCCAGATGGCTACCACCCTGGG"
def part(left, right, insert):
    return "CGTCTCA" + left + insert + right + "A" + "GAGACG"

parts = {"p1": part("AATG", "GCTT", INS1), "p2": part("GCTT", "AATG", INS2)}
result = simulate_type2s_assembly(parts, enzyme="BsmBI")
print(result["enzyme"], result["circular_product_count"], [row["length"] for row in result["products"]])
```

- **Example output:**

```text
BsmBI 1 [152]
```

## 15) `dnakit.molecular_lab.clone_into_donor_vector()`

- **Purpose:** Clone an insert into a donor vector by swapping the Type-2S stuffer.
- **Method:** Calls DnaCauldron's `swap_donor_vector_part`. The donor is circular with one site pair around its stuffer. The insert is circular, with a site pair flanking the insert, and must carry the same two overhangs as the stuffer; otherwise upstream fails and the call raises `ConfigurationError`.
- **API:** `dnakit.molecular_lab.clone_into_donor_vector(donor_sequence, insert_sequence, *, enzyme="BsmBI", verbose=True)`
- **Input:** Donor and insert DNA strings, and `enzyme` (BbsI, AarI, BsaI, BsmBI or SapI).
- **Output:** `enzyme`, `product_sequence`, `product_length` and `provenance`.
- **Example code:**

```python
from dnakit.molecular_lab import clone_into_donor_vector
INS1 = "ATGGCTAAAGGTGAAGAACTTTTCACTGGAGTTGTCCCAATTCTTGTTGAATTAGATGGTGATGTTAATGGG"
INS2 = "ATGGCTAGCGTGATCAAGCCGGAAATGCGCATCCACCTGCGTAAAAGTGTTCCAGATGGCTACCACCCTGGG"
def part(left, right, insert):
    return "CGTCTCA" + left + insert + right + "A" + "GAGACG"

donor = part("AATG", "GCTT", INS1) + "ACGT" * 175
insert = part("AATG", "GCTT", INS2)
result = clone_into_donor_vector(donor, insert, enzyme="BsmBI")
print(len(donor), result["product_length"])
```

- **Example output:**

```text
794 794
```

## 16) `dnakit.molecular_lab.assign_parts_to_backbones()`

- **Purpose:** Assign each circular Golden Gate part to the backbone whose overhangs match.
- **Method:** Calls DnaCauldron's `insert_parts_on_backbones` with BsmBI and `process_parts_with_backbone=False`. Upstream backbone detection is BsmBI only, so this wrapper accepts only BsmBI. Each part is reported as `assigned` (with its product), `no_backbone` or `error`; the per-part state is in `status_per_part`.
- **API:** `dnakit.molecular_lab.assign_parts_to_backbones(parts, backbones, *, min_backbone_length=500, verbose=True)`
- **Input:** Two mappings from name to circular DNA string: parts and backbones. Each backbone holds one BsmBI stuffer.
- **Output:** `part_count`, `assigned_count`, `assignments` (each with `part`, `status_per_part`, `already_on_backbone`, `backbone`, `product_sequence`, `product_length` and `error`), `enzyme` and `provenance`.
- **Example code:**

```python
from dnakit.molecular_lab import assign_parts_to_backbones
INS1 = "ATGGCTAAAGGTGAAGAACTTTTCACTGGAGTTGTCCCAATTCTTGTTGAATTAGATGGTGATGTTAATGGG"
INS2 = "ATGGCTAGCGTGATCAAGCCGGAAATGCGCATCCACCTGCGTAAAAGTGTTCCAGATGGCTACCACCCTGGG"
def part(left, right, insert):
    return "CGTCTCA" + left + insert + right + "A" + "GAGACG"

backbone = part("AATG", "GCTT", "ACGT" * 10) + "ACGT" * 175
parts = {"p1": part("AATG", "GCTT", INS1), "p2": part("TTAA", "CCGG", INS2)}
result = assign_parts_to_backbones(parts, {"bb": backbone})
print(result["part_count"], result["assigned_count"])
for row in result["assignments"]:
    print(row["part"], row["status_per_part"], row["backbone"])
```

- **Example output:**

```text
2 1
p1 assigned bb
p2 no_backbone None
```

## 17) `dnakit.molecular_lab.substitute_type2s_overhangs()`

- **Purpose:** Replace Type-2S overhangs on a circular vector, which moves a part within an assembly standard.
- **Method:** Calls DnaCauldron's `utils.substitute_overhangs`. Only overhangs that are sticky ends of the vector's digest are changed, and they are listed in `applied_overhangs`; `requested` lists every key given.
- **API:** `dnakit.molecular_lab.substitute_type2s_overhangs(circular_sequence, substitutions, *, enzyme="BsmBI", verbose=True)`
- **Input:** A circular DNA string; `substitutions` maps up to 8 existing 4-base overhangs to new 4-base overhangs; `enzyme` is BbsI, AarI, BsaI, BsmBI or SapI.
- **Output:** `applied_overhangs`, `requested`, `product_sequence`, `product_length`, `enzyme` and `provenance`.
- **Example code:**

```python
from dnakit.molecular_lab import substitute_type2s_overhangs
INS1 = "ATGGCTAAAGGTGAAGAACTTTTCACTGGAGTTGTCCCAATTCTTGTTGAATTAGATGGTGATGTTAATGGG"
def part(left, right, insert):
    return "CGTCTCA" + left + insert + right + "A" + "GAGACG"

donor = part("AATG", "GCTT", INS1) + "ACGT" * 175
result = substitute_type2s_overhangs(donor, {"AATG": "CCCC"}, enzyme="BsmBI")
print(result["applied_overhangs"], result["product_length"])
```

- **Example output:**

```text
['AATG'] 794
```

## 18) `dnakit.molecular_lab.record_has_backbone()`

- **Purpose:** Tell whether a circular BsmBI construct still carries a backbone besides its insert.
- **Method:** Cuts the construct with BsmBI and keeps the one site-free fragment as the insert. Calls DnaCauldron's `record_contains_backbone`; the backbone is the rest of the circle, and `has_backbone` is true when it is longer than `min_backbone_length`. Upstream ignores its enzyme argument, so the wrapper fixes BsmBI.
- **API:** `dnakit.molecular_lab.record_has_backbone(circular_sequence, *, min_backbone_length=500, verbose=True)`
- **Input:** A circular DNA string and `min_backbone_length`.
- **Output:** `construct_length`, `insert_length`, `backbone_length`, `has_backbone`, `enzyme` and `provenance`.
- **Example code:**

```python
from dnakit.molecular_lab import record_has_backbone
INS1 = "ATGGCTAAAGGTGAAGAACTTTTCACTGGAGTTGTCCCAATTCTTGTTGAATTAGATGGTGATGTTAATGGG"
def part(left, right, insert):
    return "CGTCTCA" + left + insert + right + "A" + "GAGACG"

donor = part("AATG", "GCTT", INS1) + "ACGT" * 175
result = record_has_backbone(donor)
print(result["construct_length"], result["insert_length"], result["backbone_length"], result["has_backbone"])
```

- **Example output:**

```text
794 72 722 True
```

## 19) `dnakit.molecular_lab.dna_window_overlap()`

- **Purpose:** Return the shared span of two genomic windows given as (start, end) pairs.
- **Method:** Calls DnaChisel's `windows_overlap`. Windows that only touch at one coordinate return a zero-length overlap, as DnaChisel reports it; `overlaps` is then false.
- **API:** `dnakit.molecular_lab.dna_window_overlap(window1, window2, *, verbose=True)`
- **Input:** Two (start, end) pairs of non-negative integers with start at most end.
- **Output:** `overlap` (a [start, end] pair, or None), `overlap_length`, `overlaps` and `provenance`.
- **Example code:**

```python
from dnakit.molecular_lab import dna_window_overlap
result = dna_window_overlap([0, 10], [5, 20])
print(result["overlap"], result["overlap_length"], result["overlaps"])
print(dna_window_overlap([0, 3], [5, 9])["overlap"])
```

- **Example output:**

```text
[5, 10] 5 True
None
```

## 20) `dnakit.molecular_lab.subdivide_dna_window()`

- **Purpose:** Split a (start, end) window into consecutive sub-windows of at most max_span bases.
- **Method:** Calls DnaChisel's `subdivide_window`. The pieces are contiguous and cover the whole window; an empty window gives no pieces.
- **API:** `dnakit.molecular_lab.subdivide_dna_window(window, *, max_span, verbose=True)`
- **Input:** A (start, end) pair and a positive integer `max_span`.
- **Output:** `window_count`, `windows` (each a [start, end] pair) and `provenance`.
- **Example code:**

```python
from dnakit.molecular_lab import subdivide_dna_window
result = subdivide_dna_window([0, 10], max_span=4)
print(result["window_count"], result["windows"])
```

- **Example output:**

```text
3 [[0, 4], [4, 8], [8, 10]]
```

## 21) `dnakit.molecular_lab.group_nearby_dna_segments()`

- **Purpose:** Group (start, end) segments whose start positions lie close together.
- **Method:** Sorts the segments by start and then end, and calls DnaChisel's `group_nearby_segments`. A segment joins the current group when its start is less than `max_start_gap` after the previous start and less than `max_start_spread` after the group's first start. `None` disables a test.
- **API:** `dnakit.molecular_lab.group_nearby_dna_segments(segments, *, max_start_gap=None, max_start_spread=None, verbose=True)`
- **Input:** A list of (start, end) pairs and optional positive integer limits.
- **Output:** `segment_count`, `group_count`, `groups` (lists of [start, end] pairs), the limits used and `provenance`.
- **Example code:**

```python
from dnakit.molecular_lab import group_nearby_dna_segments
segments = [[0, 5], [2, 4], [10, 12], [11, 20]]
result = group_nearby_dna_segments(segments, max_start_gap=3)
print(result["group_count"], result["groups"])
```

- **Example output:**

```text
2 [[[0, 5], [2, 4]], [[10, 12], [11, 20]]]
```

## 22) `dnakit.molecular_lab.codon_backtranslation_table()`

- **Purpose:** Return the back-translation table of a Biopython codon table as sorted codon lists.
- **Method:** Calls DnaChisel's `get_backtranslation_table`. Keys are amino-acid letters, `*` for stops, `START` for start codons and the ambiguity letters X, B, J and Z; each value is a sorted list of codons.
- **API:** `dnakit.molecular_lab.codon_backtranslation_table(table_name="Standard", *, verbose=True)`
- **Input:** `table_name`, one of Biopython's unambiguous DNA codon table names, for example `Standard` or `Vertebrate Mitochondrial`.
- **Output:** `table_name`, `entry_count`, `table` and `provenance`.
- **Example code:**

```python
from dnakit.molecular_lab import codon_backtranslation_table
result = codon_backtranslation_table()
table = result["table"]
print(result["entry_count"], table["M"], table["*"], len(table["X"]))
```

- **Example output:**

```text
26 ['ATG'] ['TAA', 'TAG', 'TGA'] 61
```

## Functions not wrapped

These rows keep a precise ledger status and a one-line reason.

- Biopython record helpers, writers and topology setters (`annotate_record`, `sequence_to_biopython_record`, `change_biopython_record_sequence`, `write_record`, `set_record_topology`, `record_is_linear`): they build or edit records rather than compute on sequences.
- Text and display helpers (`assembly2str`, `dict_to_pretty_string`, `round_all_numbers_in_dict`, `score_to_formatted_string`, the spreadsheet formatters, `gbtext_clean`): their output is not a sequence result.
- `reverse_complement_assembly` and `partially_digested_seqs_in_restriction_ligation_assembly` work on pydna's internal assembly objects; the partial-digest function returned nothing for the restriction-ligation products tried, so it has no documented example.
- `pydna.utils.express` is an upstream stub; `pydna.utils.cai` needs the `cai2` package, which is not installed in `dnakit-full` (`needs-data:cai2-package`); `parse_snapgene_history` needs a SnapGene `.dna` history file (`needs-data:snapgene-history-dna-file`).
- DnaChisel `dna_pattern_to_regexpr` maps N, R and Y to the wrong classes, so it is not wrapped. It is recorded as a duplicate of `dnakit.patterns.scan_motif`, which matches IUPAC codes correctly.
- The pandas report builders (`constraints_before_after_dataframe`, `objectives_before_after_dataframe`, `records_from_breaches_dataframe`) format results of a finished DnaChisel problem, not sequences.
