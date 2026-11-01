# Design tools: pydna, DNA Chisel and DnaCauldron

This group of 24 functions in `dnakit.molecular_lab` fills gaps left by the first set of adapters over three construct-design packages: [pydna](https://github.com/pydna-group/pydna) (gel prediction and drawing, codon adaptation index, SnapGene files, NCBI BLAST and the NEB Tm service), [DnaChisel](https://github.com/Edinburgh-Genome-Foundry/DnaChisel) (sequence composition, pattern annotation, constraint-breach tables and plots, optimization reports, record loading, BLAST and Bowtie searches) and [DnaCauldron](https://github.com/Edinburgh-Genome-Foundry/DnaCauldron) (assembly-plan simulation, a leveled plan graph and a PDF summary). The functions are imported as `from dnakit.molecular_lab import ...`. None of them designs or predicts a wet-lab result beyond what the upstream routine computes.

The upstream packages are optional dependencies and are imported only when a function is called; a missing package raises `BackendUnavailableError` with the install command. Licences, read from the installed package metadata: pydna is BSD, DnaChisel is MIT and DnaCauldron is MIT. No GPL package is involved. Three functions call external programs through DNAKit's bounded command runner ([BLAST+](https://blast.ncbi.nlm.nih.gov/doc/blast-help/downloadblastdata.html) `blastn`, and [Bowtie 1](https://github.com/BenLangmead/bowtie) `bowtie` and `bowtie-build`), and two wrap network services (NCBI BLAST and the NEB Tm calculator). The report and plot functions additionally need `pdf_reports` and `sequenticon`, `geneblocks`, or WeasyPrint, as each section states. Every function that writes a file refuses to overwrite an existing path.

## Example data

Two examples read declared input files: a one-record FASTA file with a `GGTCTC` site at positions 3 to 9, and a GenBank record whose constraint is annotated in a feature label. All other examples pass their inputs directly or write their files into a temporary directory.

**Example file:** `part.fa`

```text
>part1
ATGGGTCTCATGCCCAAAGCTTTAA
```

**Example file:** `annotated.gb`

```text
LOCUS       annotated                 25 bp    DNA     linear   UNK 01-JAN-1980
DEFINITION  .
FEATURES             Location/Qualifiers
     misc_feature    4..9
                     /label="@AvoidPattern(GGTCTC)"
ORIGIN
        1 atgggtctca tgcccaaagc tttaa
//
```

## 1) `dnakit.molecular_lab.pydna_blast_best_nt_hit()`

- **Purpose:** Find the best NCBI nt hit for a DNA sequence with BLAST and return the matched region.
- **Method:** Validate the sequence, then call pydna's `sequence_picker.genbank_accession`, which submits a megablast search against NCBI's nt database and returns the best hit region as a record. Despite its name it looks up no accession. `name` is only echoed back as `query_name`; it is not sent to BLAST.
- **API:** `pydna_blast_best_nt_hit(sequence: 'str', *, name: 'str' = 'query', verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the optional package pydna (extra `dnakit[bio-assembly]`) and internet access to the NCBI BLAST service (no key). It cannot run offline, so the example shows the input check that runs before any network call. Usage with the dependency installed is shown after the example.
- **Input:** `sequence` is a string of A, C, G and T (whitespace is removed, case is ignored) of at most 20000 bases. `name` is a non-empty string of at most 128 characters.
- **Output:** `service`, `source` (the BLAST URL) and `result` with `query_name`, `sequence` (the matched region, upper case), `length`, `circular` and `description`.
- **Example code:**

```python
from dnakit.exceptions import ConfigurationError
from dnakit.molecular_lab import pydna_blast_best_nt_hit

try:
    pydna_blast_best_nt_hit("ATGCNNGTAC", verbose=False)
except ConfigurationError as error:
    print(type(error).__name__, "only A, C, G and T" in str(error))
```

- **Example output:**

```text
ConfigurationError True
```

Usage with the dependency installed:

```python
from dnakit.molecular_lab import pydna_blast_best_nt_hit

hit = pydna_blast_best_nt_hit("ATGGCTAAAGCT" * 6, name="insert", verbose=False)
print(hit["result"]["length"], hit["result"]["description"])
```

## 2) `dnakit.molecular_lab.pydna_codon_adaptation_index()`

- **Purpose:** Compute the codon adaptation index (CAI) of a coding sequence.
- **Method:** With explicit `weights` the index is computed in DNAKit itself: the geometric mean of the weights of the scored codons (Sharp and Li, 1987), rounded to 3 decimals as pydna does. ATG, TGG and the stop codons TAA, TAG and TGA are not scored. Without `weights`, the function calls pydna's `utils.cai` with the organism table, which needs the cai2 package.
- **API:** `pydna_codon_adaptation_index(seq: 'str', *, weights: 'dict[str, float] | None' = None, organism: 'str' = 'sce', verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** With `weights` it works after plain pip install (it imports no optional package). Without `weights` it needs pydna (extra `dnakit[bio-assembly]`) and cai2 (`pip install "pydna[express]"`).
- **Input:** `seq` is a string of A, C, G and T, a whole number of codons, at most 300000 bases. `weights` maps three-letter codons to positive finite numbers; every scored codon of `seq` must have a weight, otherwise `ConfigurationError` lists the missing ones. The docstring says weights lie in (0, 1], but the code only requires them to be positive. `organism` (default `sce`) must be a key of pydna's codon table and is used only without `weights`.
- **Output:** `status`, `operation`, `method`, `cai`, `codon_count`, `scored_codon_count` (explicit weights only) and `provenance`.
- **Example code:**

```python
from dnakit.molecular_lab import pydna_codon_adaptation_index

weights = {"GCT": 1.0, "GCC": 0.5, "AAA": 0.8, "AAG": 0.4}
result = pydna_codon_adaptation_index("ATGGCTGCCAAAAAGTAA", weights=weights, verbose=False)
print(result["cai"], result["codon_count"], result["scored_codon_count"])
```

- **Example output:**

```text
0.632 6 4
```

## 3) `dnakit.molecular_lab.pydna_draw_gel()`

- **Purpose:** Draw an agarose gel image from fragment lengths and write it as a new PNG.
- **Method:** Build pydna's gel size model (a cubic spline through the migration of the chosen ladder), turn each fragment length into a `FakeSeq`, call `pydna.gel.gel` and save the Pillow image as PNG. Every band length must lie inside the ladder range.
- **API:** `pydna_draw_gel(output_path: 'str', *, lanes: 'list[list[int]] | None' = None, ladder: 'str' = 'GeneRuler_1kb_plus', gel_length: 'int' = 600, margin: 'int' = 50, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the optional package pydna (extra `dnakit[bio-assembly]`) with its gel extra, which provides Pillow and SciPy (`pip install "pydna[gel]"`).
- **Input:** `output_path` ends in `.png` and does not exist. `lanes` is a list of 1 to 20 lanes, each a list of 1 to 50 fragment lengths in bp; with `None` the ladder itself is drawn as one lane. `ladder` is `GeneRuler_1kb_plus` or `HI_LO_DNA_MARKER`. `gel_length` is 100 to 5000 and `margin` is 0 to `gel_length - 10`.
- **Output:** `path` and `data` with `ladder`, `lanes`, `lane_count`, `band_count`, `gel_length`, `margin`, `image_width`, `image_height` and `provenance`.
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.molecular_lab import pydna_draw_gel

with tempfile.TemporaryDirectory() as folder:
    result = pydna_draw_gel(
        str(Path(folder) / "gel.png"),
        lanes=[[500, 1000], [2000, 3000, 5000]],
        verbose=False,
    )
    data = result["data"]
    print(data["lane_count"], data["band_count"], data["image_width"], data["image_height"])
    print(Path(result["path"]).read_bytes()[:4])
```

- **Example output:**

```text
2 5 180 600
b'\x89PNG'
```

## 4) `dnakit.molecular_lab.pydna_neb_melting_temperature()`

- **Purpose:** Return the melting temperature (degrees Celsius) of a primer from the NEB Tm calculator.
- **Method:** Validate the primer, concentration and product code, then call pydna's `tm.tm_neb`, which posts the primer to the NEB service at `https://tmapi.neb.com/tm`.
- **API:** `pydna_neb_melting_temperature(primer: 'str', *, conc: 'float' = 0.5, prodcode: 'str' = 'q5-0', verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the optional package pydna (extra `dnakit[bio-assembly]`) and internet access to the NEB Tm API (no key). It cannot run offline, so the example shows the input check that runs before any network call. Usage with the dependency installed is shown after the example.
- **Input:** `primer` is a string of A, C, G and T of at most 200 bases. `conc` is the primer concentration in micromolar, 0.001 to 100. `prodcode` is an NEB product code of 1 to 64 letters, digits, `-` or `_` (for example `q5-0` for Q5 polymerase).
- **Output:** `service`, `source` and `result` with `primer`, `conc_uM`, `prodcode` and `tm_celsius`.
- **Example code:**

```python
from dnakit.exceptions import ConfigurationError
from dnakit.molecular_lab import pydna_neb_melting_temperature

try:
    pydna_neb_melting_temperature("ATGCGTACGTTAGC", prodcode="q5 0", verbose=False)
except ConfigurationError as error:
    print(type(error).__name__, "prodcode must be" in str(error))
```

- **Example output:**

```text
ConfigurationError True
```

Usage with the dependency installed:

```python
from dnakit.molecular_lab import pydna_neb_melting_temperature

result = pydna_neb_melting_temperature("ATGCGTACGTTAGC", conc=0.5, prodcode="q5-0", verbose=False)
print(result["result"]["tm_celsius"])
```

## 5) `dnakit.molecular_lab.pydna_parse_snapgene_history()`

- **Purpose:** Parse a local SnapGene `.dna` file with pydna, resolving its cloning history when it has one.
- **Method:** Call pydna's `snapgene_history_parser.parse_snapgene_history` on the file, summarize the root record and its features, and ask `sgffp.SgffReader` whether the file carries a history. A file without a history returns the root sequence only.
- **API:** `pydna_parse_snapgene_history(path: 'str', *, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the optional package pydna (extra `dnakit[bio-assembly]`), which depends on sgffp.
- **Input:** `path` is a local path ending in `.dna`; URLs are rejected.
- **Output:** `name`, `length`, `circular`, `sequence` (upper case), `feature_count`, `features` (each with `type`, `start`, `end`, `strand`, `label`), `has_history`, `source_type` (the class name of the record's source, or `None`) and `provenance`.
- **Example code:**

```python
import tempfile
from pathlib import Path

import sgffp

from dnakit.molecular_lab import pydna_parse_snapgene_history

sequence = "ATGCGTACGTTAGCATGCGTACGTTAGC"
with tempfile.TemporaryDirectory() as folder:
    path = Path(folder) / "construct.dna"
    obj = sgffp.SgffObject.new(sequence=sequence, topology="circular")
    path.write_bytes(sgffp.SgffWriter.to_bytes(obj))
    result = pydna_parse_snapgene_history(str(path), verbose=False)
print(result["name"], result["length"], result["circular"])
print(result["has_history"], result["feature_count"], result["sequence"] == sequence)
```

- **Example output:**

```text
construct 28 True
False 0 True
```

## 6) `dnakit.molecular_lab.pydna_predict_gel_relative_front()`

- **Purpose:** Predict the relative migration front (Rf) of DNA fragments on an agarose gel.
- **Method:** Build pydna's gel size model, a natural cubic spline through the Rf values of a built-in ladder (`pydna.gel.interpolator`), and evaluate it at each fragment length. Lengths outside the ladder range have no prediction (`None`). No figure is drawn; use `pydna_draw_gel` for an image.
- **API:** `pydna_predict_gel_relative_front(lengths: 'list[int]', *, ladder: 'str' = 'GeneRuler_1kb_plus', verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the optional package pydna (extra `dnakit[bio-assembly]`) with SciPy (`pip install "pydna[gel]"`).
- **Input:** `lengths` is a list of 1 to 500 integers between 1 and 10000000 (bp). `ladder` is `GeneRuler_1kb_plus` or `HI_LO_DNA_MARKER`.
- **Output:** `ladder`, `ladder_sizes`, `ladder_relative_fronts`, `lengths`, `relative_front` (one float or `None` per length) and `provenance`.
- **Example code:**

```python
from dnakit.molecular_lab import pydna_predict_gel_relative_front

result = pydna_predict_gel_relative_front([500, 1000, 3000, 50000], verbose=False)
print(result["ladder"], len(result["ladder_sizes"]))
print([None if rf is None else round(rf, 3) for rf in result["relative_front"]])
```

- **Example output:**

```text
GeneRuler_1kb_plus 15
[0.755, 0.562, 0.261, None]
```

## 7) `dnakit.molecular_lab.annotate_pattern_occurrences_features()`

- **Purpose:** Find every occurrence of an IUPAC DNA pattern on both strands and annotate them.
- **Method:** Turn each IUPAC code of the pattern into the full set of bases it stands for and scan the sequence with a look-ahead regular expression, so overlapping hits are kept. The reverse complement of the pattern is scanned for strand -1 unless the pattern is palindromic (then only the forward strand is searched, as DnaChisel does). Optionally write the features to a new GenBank file with Biopython. The function says that DnaChisel's own pattern regex maps N, R and Y to the wrong sets, so its annotations differ for patterns containing those codes.
- **API:** `annotate_pattern_occurrences_features(sequence: 'str', pattern: 'str', *, feature_type: 'str' = 'misc_feature', prefix: 'str' = '!', output_path: 'str | Path | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Works after plain pip install. Writing `output_path` needs Biopython (extra `dnakit[bio-formats]`).
- **Input:** `sequence` and `pattern` are IUPAC DNA strings (the sequence may be empty; the pattern is at most 1000 bases). `feature_type` is a GenBank feature key; `prefix` (at most 20 characters) is prepended to the pattern to form the label. `output_path` must not exist and must end in `.gb`, `.gbk` or `.genbank`.
- **Output:** `sequence_length`, `pattern`, `palindromic`, `occurrence_count`, `features` (each with `type`, `start` inclusive, `end` exclusive, `strand` +1 or -1, `label`), `output_path` and `provenance`.
- **Example code:**

```python
from dnakit.molecular_lab import annotate_pattern_occurrences_features

result = annotate_pattern_occurrences_features("AAGGTCTCTTGAGACCAA", "GGTCTC", verbose=False)
print(result["palindromic"], result["occurrence_count"])
for feature in result["features"]:
    print(feature["start"], feature["end"], feature["strand"], feature["label"])
```

- **Example output:**

```text
False 2
2 8 1 !GGTCTC
10 16 -1 !GGTCTC
```

## 8) `dnakit.molecular_lab.blast_sequence_hits()`

- **Purpose:** BLAST one sequence against a BLAST database or a subject FASTA file with `blastn`.
- **Method:** Write the query to a temporary FASTA file, run `blastn -outfmt 5` with DnaChisel's defaults (word size 4, 80 % identity, up to 1000 alignments, dust off) through the bounded runner (no shell), and parse the XML with Biopython's `Bio.Blast.NCBIXML` into alignments and high-scoring pairs. The `blastn` version is read first with `-version`.
- **API:** `blast_sequence_hits(sequence: 'str', *, blast_db: 'str | None' = None, subject_fasta: 'str | Path | None' = None, word_size: 'int' = 4, perc_identity: 'float' = 80.0, num_alignments: 'int' = 1000, ungapped: 'bool' = False, num_threads: 'int' = 3, culling_limit: 'int | None' = None, e_value: 'float | None' = None, task: 'str | None' = None, dust: 'str' = 'no', blastn_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs Biopython (extra `dnakit[bio-formats]`) and the external program `blastn` (NCBI BLAST+) on PATH, or passed as `blastn_path`. `blastn` is not installed in the environment used to build this page, so the example points `blastn_path` at a file that does not exist and shows the resulting error; this gives the same output on every machine. Usage with the program installed is shown after the example.
- **Input:** `sequence` is an IUPAC DNA string. Exactly one of `blast_db` (a database prefix made with `makeblastdb`, not starting with `-`) or `subject_fasta` (an existing file) is required. `word_size` is 4 to 28, `perc_identity` 0 to 100, `num_alignments` 1 to 10000, `num_threads` 1 to 64, `task` is one of `blastn`, `blastn-short`, `dc-megablast`, `megablast`, `dust` is `no`, `yes` or three integers such as `20 64 1`.
- **Output:** `status`, `operation`, `backend`, `version`, `blast_db`, `subject_fasta`, `queries` (each with `query_id`, `query_length` and `alignments`; each alignment has `title`, `length` and `hsps` with `identities`, `align_length`, `gaps`, `query_start`, `query_end`, `subject_start`, `subject_end`, `score`, `evalue`), `hit_count` and `provenance`.
- **Example code:**

```python
from dnakit.exceptions import BackendUnavailableError
from dnakit.molecular_lab import blast_sequence_hits

try:
    blast_sequence_hits(
        "ATGGGTCTCATGCCCAAA",
        blast_db="nt_database",
        blastn_path="/no/such/dir/blastn",
        verbose=False,
    )
except BackendUnavailableError as error:
    print(type(error).__name__, "missing or not executable" in str(error))
```

- **Example output:**

```text
BackendUnavailableError True
```

Usage with the dependency installed:

```python
from dnakit.molecular_lab import blast_sequence_hits

result = blast_sequence_hits("ATGGGTCTCATGCCCAAA", subject_fasta="genes.fa", verbose=False)
print(result["hit_count"])
```

## 9) `dnakit.molecular_lab.constraints_breaches_table()`

- **Purpose:** Return a table of constraint breaches for several sequences, one row per sequence.
- **Method:** Build the DnaChisel constraints from the allow-list used by `evaluate_dnachisel_specs` and call DnaChisel's `constraints_breaches_dataframe`. Each column is one constraint, labelled in DnaChisel's breach form; each cell lists the breached locations as `start-end(strand)` separated by commas, and is empty when the constraint passes.
- **API:** `constraints_breaches_table(sequences: 'object', *, constraints: 'Sequence[Mapping[str, Any]]', display_constraints_locations: 'bool' = False, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the optional package DnaChisel (extra `dnakit[bio-dnachisel]`).
- **Input:** `sequences` is a mapping `{name: sequence}` or a list of `{"name", "sequence"}` mappings (1 to 1000 unique names, IUPAC DNA, at most 2000000 bases in total). `constraints` is a non-empty list of mappings with a `name` from the allow-list (`AvoidPattern`, `AvoidChanges`, `AvoidHairpins`, `AvoidStopCodons`, `CodonOptimize`, `EnforceChanges`, `EnforceGCContent`, `EnforceMeltingTemperature`, `EnforceTranslation`, `MaximizeCAI`) and keyword arguments. `display_constraints_locations` is passed to DnaChisel and changes the column labels.
- **Output:** `sequence_count`, `constraint_labels`, `rows` (one dictionary per sequence with its `sequence` name and one key per label) and `provenance`.
- **Example code:**

```python
from dnakit.molecular_lab import constraints_breaches_table

constraints = [
    {"name": "AvoidPattern", "pattern": "GGTCTC"},
    {"name": "EnforceGCContent", "mini": 0.3, "maxi": 0.7, "window": 18},
]
result = constraints_breaches_table(
    {"a": "ATGGGTCTCATGCCCAAA", "b": "ATGGAGACCAAGCCCAAA", "c": "A" * 22},
    constraints=constraints,
    verbose=False,
)
print(result["constraint_labels"])
for row in result["rows"]:
    print(row)
```

- **Example output:**

```text
['GGTCTC', 'GC outside 30-70%/18bp']
{'sequence': 'a', 'GGTCTC': '3-9(+)', 'GC outside 30-70%/18bp': ''}
{'sequence': 'b', 'GGTCTC': '3-9(-)', 'GC outside 30-70%/18bp': ''}
{'sequence': 'c', 'GGTCTC': '', 'GC outside 30-70%/18bp': '0-22'}
```

## 10) `dnakit.molecular_lab.create_bowtie_index()`

- **Purpose:** Build a Bowtie 1 index from sequences; sequences are named by their position (0, 1, ...).
- **Method:** Create the new directory `index_dir`, write the sequences to `sequences.fa` in it and run `bowtie-build -f sequences.fa bowtie --quiet` through the bounded runner (no shell). The index files are `bowtie.*.ebwt`; the returned `index_prefix` is what `find_bowtie_matches` takes.
- **API:** `create_bowtie_index(sequences: 'Sequence[str]', *, index_dir: 'str | Path', bowtie_build_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the external program `bowtie-build` (Bowtie 1) on PATH, or passed as `bowtie_build_path`. It is not installed in the environment used to build this page, so the example points `bowtie_build_path` at a file that does not exist and shows the resulting error. Usage with the program installed is shown after the example.
- **Input:** `sequences` is a non-empty list of at most 1000 IUPAC DNA strings (at most 10000000 bases in total). `index_dir` must not exist yet.
- **Output:** `status`, `operation`, `backend`, `version`, `index_prefix`, `index_files` (relative names), `sequence_count` and `provenance`.
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.exceptions import BackendUnavailableError
from dnakit.molecular_lab import create_bowtie_index

with tempfile.TemporaryDirectory() as folder:
    try:
        create_bowtie_index(
            ["ATGGGTCTCATGCCCAAAGCTTTAA"],
            index_dir=Path(folder) / "index",
            bowtie_build_path="/no/such/dir/bowtie-build",
            verbose=False,
        )
    except BackendUnavailableError as error:
        print(type(error).__name__, "missing or not executable" in str(error))
    print((Path(folder) / "index").exists())
```

- **Example output:**

```text
BackendUnavailableError True
False
```

Usage with the dependency installed:

```python
from dnakit.molecular_lab import create_bowtie_index

index = create_bowtie_index(["ATGGGTCTCATGCCCAAAGCTTTAA"], index_dir="bowtie_index", verbose=False)
print(index["index_prefix"], index["sequence_count"])
```

## 11) `dnakit.molecular_lab.find_bowtie_matches()`

- **Purpose:** Find, for every k-mer of a sequence, the best Bowtie 1 hit in an index (DnaChisel's matcher).
- **Method:** Cut the sequence into k-mers of `match_length`, write them to a temporary FASTA file and run `bowtie --best -k 1 -v <max_mismatches> <index> -f ...` through the bounded runner. Each reported alignment becomes `{"start", "end", "mismatches"}`, where `start` and `end` give the k-mer position in `sequence` and `mismatches` counts the mismatches in Bowtie's descriptor. The function checks that `<index_prefix>.1.ebwt` exists before it looks for the program.
- **API:** `find_bowtie_matches(sequence: 'str', *, index_prefix: 'str | Path', match_length: 'int', max_mismatches: 'int' = 0, bowtie_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the external program `bowtie` (Bowtie 1) on PATH, or passed as `bowtie_path`, and an index built by `create_bowtie_index`. It is not installed in the environment used to build this page, so the example uses a placeholder index file and a `bowtie_path` that does not exist, and shows the resulting error. Usage with the program installed is shown after the example.
- **Input:** `sequence` is an IUPAC DNA string. `index_prefix` is the prefix of an existing index. `match_length` is 1 to the length of `sequence`. `max_mismatches` is 0 to 3.
- **Output:** `status`, `operation`, `backend`, `version`, `kmer_length`, `match_count`, `matches` (sorted by `start`) and `provenance`.
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.exceptions import BackendUnavailableError
from dnakit.molecular_lab import find_bowtie_matches

with tempfile.TemporaryDirectory() as folder:
    prefix = Path(folder) / "bowtie"
    Path(f"{prefix}.1.ebwt").write_bytes(b"")  # placeholder; only its existence is checked first
    try:
        find_bowtie_matches(
            "ATGGGTCTCATGCCCAAA",
            index_prefix=prefix,
            match_length=8,
            bowtie_path="/no/such/dir/bowtie",
            verbose=False,
        )
    except BackendUnavailableError as error:
        print(type(error).__name__, "missing or not executable" in str(error))
```

- **Example output:**

```text
BackendUnavailableError True
```

Usage with the dependency installed:

```python
from dnakit.molecular_lab import create_bowtie_index, find_bowtie_matches

index = create_bowtie_index(["ATGGGTCTCATGCCCAAAGCTTTAA"], index_dir="bowtie_index", verbose=False)
found = find_bowtie_matches(
    "ATGGGTCTCATGCCC", index_prefix=index["index_prefix"], match_length=8, verbose=False
)
print(found["match_count"], found["matches"][:2])
```

## 12) `dnakit.molecular_lab.gc_fraction_or_windows()`

- **Purpose:** Return the global GC fraction, or the GC fraction of every sliding window.
- **Method:** Upper-case the sequence and call DnaChisel's `gc_content`. The global value is the count of G and C divided by the sequence length, where every symbol including N counts in the denominator. With `window_size` the result has one value per window start `i` (`0 <= i <= len - window_size`), covering `[i, i + window_size)`.
- **API:** `gc_fraction_or_windows(sequence: 'str', *, window_size: 'int | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the optional package DnaChisel (extra `dnakit[bio-dnachisel]`).
- **Input:** `sequence` is an IUPAC DNA string; it may be empty only when `window_size` is given. `window_size` is an integer from 1 up to 200000; a window longer than the sequence gives no windows.
- **Output:** `sequence_length`, `window_size`, `gc_fraction` (`None` in window mode), `window_count`, `window_gc_fractions` and `provenance`.
- **Example code:**

```python
from dnakit.molecular_lab import gc_fraction_or_windows

whole = gc_fraction_or_windows("ATGCGCNNAT", verbose=False)
print(whole["sequence_length"], whole["gc_fraction"])
windows = gc_fraction_or_windows("ATGCGCAT", window_size=4, verbose=False)
print(windows["window_count"], windows["window_gc_fractions"])
```

- **Example output:**

```text
10 0.4
5 [0.5, 0.75, 1.0, 0.75, 0.5]
```

## 13) `dnakit.molecular_lab.generate_random_dna_with_composition()`

- **Purpose:** Generate a random DNA sequence with a given GC share or per-base probabilities.
- **Method:** Call DnaChisel's `random_dna_sequence`. `gc_share` is the total G+C probability, split evenly between G and C and between A and T. `probas` maps bases to probabilities and must sum to 1. DnaChisel draws with numpy's global generator, so `seed` seeds that generator (a side effect of the original routine).
- **API:** `generate_random_dna_with_composition(length: 'int', *, gc_share: 'float | None' = None, probas: 'Mapping[str, float] | None' = None, seed: 'int | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the optional package DnaChisel (extra `dnakit[bio-dnachisel]`).
- **Input:** `length` is 0 to 200000. `gc_share` is 0 to 1. `probas` has distinct keys among A, C, G and T. Giving both `gc_share` and `probas` is an error; giving neither makes the four bases equiprobable. `seed` is 0 to 4294967295.
- **Output:** `sequence`, `length`, `composition_source` (`gc_share`, `probas` or `uniform`), `gc_share`, `probas`, `gc_fraction_observed` (`None` for length 0), `seed` and `provenance`.
- **Example code:**

```python
from dnakit.molecular_lab import generate_random_dna_with_composition

first = generate_random_dna_with_composition(20, gc_share=0.5, seed=1, verbose=False)
second = generate_random_dna_with_composition(20, gc_share=0.5, seed=1, verbose=False)
print(first["sequence"], first["composition_source"], first["gc_fraction_observed"])
print(first["sequence"] == second["sequence"])
```

- **Example output:**

```text
CAGCGGGCCACAGTGACAGG gc_share 0.7
True
```

## 14) `dnakit.molecular_lab.list_common_restriction_enzymes_filtered()`

- **Purpose:** List restriction enzymes by site length, optimum temperature, supplier count and unlike sites.
- **Method:** Call DnaChisel's `list_common_enzymes` with its full parameter set (DNAKit's `list_common_restriction_enzymes` keeps only a single site length and the supplier count) and sort the names. `site_unlike` names enzymes whose recognition sites, in all IUPAC variants, must not occur among the returned enzymes' sites. The table is DnaChisel's bundled enzyme table (Biopython Restriction), so the exact lists depend on the installed versions.
- **API:** `list_common_restriction_enzymes_filtered(*, site_length: 'Sequence[int]' = (6,), opt_temp: 'Sequence[int]' = (37,), min_suppliers: 'int' = 1, site_unlike: 'Sequence[str]' = (), verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the optional package DnaChisel (extra `dnakit[bio-dnachisel]`).
- **Input:** `site_length` is 1 to 20 integers from 4 to 12. `opt_temp` is 1 to 20 integers from 0 to 120 (degrees Celsius). `min_suppliers` is 0 to 20. `site_unlike` is a list of enzyme names.
- **Output:** `site_length`, `opt_temp`, `min_suppliers`, `site_unlike`, `enzyme_count`, `enzymes` (sorted names) and `provenance`.
- **Example code:**

```python
from dnakit.molecular_lab import list_common_restriction_enzymes_filtered

six = list_common_restriction_enzymes_filtered(verbose=False)
eight = list_common_restriction_enzymes_filtered(site_length=[8], verbose=False)
print("EcoRI" in six["enzymes"], "EcoRI" in eight["enzymes"], "NotI" in eight["enzymes"])
print(six["enzymes"] == sorted(six["enzymes"]), six["enzyme_count"] == len(six["enzymes"]))
```

- **Example output:**

```text
True False True
True True
```

## 15) `dnakit.molecular_lab.load_dna_record()`

- **Purpose:** Read one sequence record from FASTA, GenBank, EMBL or SnapGene (`.dna`) and return it.
- **Method:** With `file_format="auto"` choose the reader by extension: `.gb`, `.gbk` and `.genbank` are GenBank, `.fa`, `.fasta` and `.fna` are FASTA, and `.dna` is SnapGene. GenBank, EMBL and FASTA go through `Bio.SeqIO.read`, so a file with more than one record is an error; SnapGene files go through `snapgene_reader`. The sequence keeps the case used in the file.
- **API:** `load_dna_record(source: 'str | Path', *, file_format: 'str' = 'auto', name: 'str | None' = None, linear: 'bool' = True, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs Biopython (extra `dnakit[bio-formats]`) for FASTA, GenBank and EMBL, and snapgene_reader (extra `dnakit[bio-snapgene]`) for `.dna` files.
- **Input:** `source` is an existing file of at most 50 MB. `file_format` is `auto`, `genbank`, `embl`, `fasta` or `snapgene`. `name`, if given, replaces the record id and name. `linear` sets only the reported topology; the files are not inspected for it.
- **Output:** `format`, `id`, `name`, `description`, `sequence`, `length`, `topology`, `feature_count`, `source` (the path as given), `status`, `operation` and `provenance`.
- **Example code:**

```python
from dnakit.molecular_lab import load_dna_record

record = load_dna_record("part.fa", verbose=False)
print(record["format"], record["id"], record["length"], record["topology"])
print(record["sequence"])
```

- **Example output:**

```text
fasta part1 25 linear
ATGGGTCTCATGCCCAAAGCTTTAA
```

## 16) `dnakit.molecular_lab.plot_constraint_breaches_png()`

- **Purpose:** Draw one sequence with its constraint breaches highlighted and save it as a PNG.
- **Method:** Build the constraint-breach table for the single sequence, convert it with DnaChisel's `records_from_breaches_dataframe` and draw it with `plot_breaches_record` (matplotlib and DNA Features Viewer). The breached locations are parsed from the table and returned next to the image path.
- **API:** `plot_constraint_breaches_png(sequence: 'str', *, output_path: 'str | Path', constraints: 'Sequence[Mapping[str, Any]]', name: 'str' = 'sequence', figure_width: 'float' = 10.0, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the optional packages DnaChisel (extra `dnakit[bio-dnachisel]`), matplotlib (extra `dnakit[bio-plots]`) and DNA Features Viewer (extra `dnakit[bio-dna-features-viewer]`).
- **Input:** `sequence` is an IUPAC DNA string. `output_path` must not exist and must end in `.png`. `constraints` is a non-empty list of specification mappings as in section 9. `figure_width` is 2 to 40.
- **Output:** `image_path`, `sequence_length`, `breaches` (for each constraint label a list of `[start, end, strand]`, with strand `1`, `-1` or `None`), `constraint_labels`, `status`, `operation` and `provenance`.
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.molecular_lab import plot_constraint_breaches_png

constraints = [
    {"name": "AvoidPattern", "pattern": "GGTCTC"},
    {"name": "EnforceGCContent", "mini": 0.3, "maxi": 0.7, "window": 18},
]
with tempfile.TemporaryDirectory() as folder:
    result = plot_constraint_breaches_png(
        "ATGGGTCTCATGCCCAAAGAGACCAA",
        output_path=Path(folder) / "breaches.png",
        constraints=constraints,
        verbose=False,
    )
    print(Path(result["image_path"]).read_bytes()[:4])
print(result["breaches"])
```

- **Example output:**

```text
b'\x89PNG'
{'GGTCTC': [[3, 9, 1], [18, 24, -1]], 'GC outside 30-70%/18bp': []}
```

## 17) `dnakit.molecular_lab.plot_optimization_changes_png()`

- **Purpose:** Optimize a sequence and save a PNG of the changed regions (DnaChisel's `plot_optimization_changes`).
- **Method:** Build the problem, resolve the constraints (and optimize the objectives when there are any), then call DnaChisel's `plot_optimization_changes` and save the figure. The changed segments are computed by comparing the sequence before and after.
- **API:** `plot_optimization_changes_png(sequence: 'str', *, output_path: 'str | Path', constraints: 'Sequence[Mapping[str, Any]]' = (), objectives: 'Sequence[Mapping[str, Any]]' = (), circular: 'bool' = False, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the optional packages DnaChisel (extra `dnakit[bio-dnachisel]`), geneblocks (extra `dnakit[bio-geneblocks]`, which also installs DNA Features Viewer) and matplotlib (extra `dnakit[bio-plots]`). geneblocks is not installed in the environment used to build this page, so the example shows the resulting `BackendUnavailableError`. Usage with the dependency installed is shown after the example.
- **Input:** `sequence` is an IUPAC DNA string. `output_path` must not exist and must end in `.png`. `constraints` and `objectives` are lists of specification mappings as in section 9; at least one of the two must be non-empty. `circular=True` treats the sequence as circular. If the constraints cannot be satisfied, `ConfigurationError` is raised.
- **Output:** `image_path`, `sequence_length`, `changed_segments` (half-open `[start, end)`), `edits`, `all_constraints_pass`, `status`, `operation` and `provenance`.
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.exceptions import BackendUnavailableError
from dnakit.molecular_lab import plot_optimization_changes_png

with tempfile.TemporaryDirectory() as folder:
    try:
        plot_optimization_changes_png(
            "ATGGGTCTCATGCCCAAA",
            output_path=Path(folder) / "changes.png",
            constraints=[{"name": "AvoidPattern", "pattern": "GGTCTC"}],
            verbose=False,
        )
    except BackendUnavailableError as error:
        print(type(error).__name__, "geneblocks could not be imported" in str(error))
```

- **Example output:**

```text
BackendUnavailableError True
```

Usage with the dependency installed:

```python
from dnakit.molecular_lab import plot_optimization_changes_png

result = plot_optimization_changes_png(
    "ATGGGTCTCATGCCCAAA",
    output_path="changes.png",
    constraints=[{"name": "AvoidPattern", "pattern": "GGTCTC"}],
    verbose=False,
)
print(result["changed_segments"], result["edits"])
```

## 18) `dnakit.molecular_lab.write_constraint_breaches_pdf()`

- **Purpose:** Write one breach page per sequence to a new PDF and return the breach table as well.
- **Method:** Build the breach table for all sequences, convert it with `records_from_breaches_dataframe` and write the pages with DnaChisel's `breaches_records_to_pdf` (matplotlib and DNA Features Viewer).
- **API:** `write_constraint_breaches_pdf(sequences: 'object', *, output_path: 'str | Path', constraints: 'Sequence[Mapping[str, Any]]', figure_width: 'float' = 10.0, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the optional packages DnaChisel (extra `dnakit[bio-dnachisel]`), matplotlib (extra `dnakit[bio-plots]`) and DNA Features Viewer (extra `dnakit[bio-dna-features-viewer]`).
- **Input:** `sequences` and `constraints` as in section 9. `output_path` must not exist and must end in `.pdf`. `figure_width` is 2 to 40.
- **Output:** `pdf_path`, `page_count`, `sequence_count`, `constraint_labels`, `rows` (the breach table), `status`, `operation` and `provenance`.
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.molecular_lab import write_constraint_breaches_pdf

constraints = [{"name": "AvoidPattern", "pattern": "GGTCTC"}]
with tempfile.TemporaryDirectory() as folder:
    result = write_constraint_breaches_pdf(
        {"a": "ATGGGTCTCATGCCCAAA", "b": "ATGGAGACCAAGCCCAAA"},
        output_path=Path(folder) / "breaches.pdf",
        constraints=constraints,
        verbose=False,
    )
    print(Path(result["pdf_path"]).read_bytes()[:4])
print(result["page_count"], result["rows"])
```

- **Example output:**

```text
b'%PDF'
2 [{'sequence': 'a', 'GGTCTC': '3-9(+)'}, {'sequence': 'b', 'GGTCTC': '3-9(-)'}]
```

## 19) `dnakit.molecular_lab.write_dnachisel_no_solution_report()`

- **Purpose:** Report where constraints cannot be satisfied, or report that they can.
- **Method:** Build the problem and run `resolve_constraints()` (and `optimize()` when objectives are given). If DnaChisel raises `NoSolutionError` with a location, call DnaChisel's `write_no_solution_report` to write `plots.pdf`, the local breach GenBank file and the logs into the new `output_dir`, and return the breach location. If the error has no location, the status is `no_solution_without_location` and nothing is written. If the constraints are satisfiable the status is `solved` and nothing is written.
- **API:** `write_dnachisel_no_solution_report(sequence: 'str', *, output_dir: 'str | Path', constraints: 'Sequence[Mapping[str, Any]]' = (), objectives: 'Sequence[Mapping[str, Any]]' = (), circular: 'bool' = False, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the optional packages DnaChisel (extra `dnakit[bio-dnachisel]`), matplotlib (extra `dnakit[bio-plots]`) and DNA Features Viewer (extra `dnakit[bio-dna-features-viewer]`). It does not need the report packages `pdf_reports` and `sequenticon`.
- **Input:** `sequence` is an IUPAC DNA string. `output_dir` must not exist (its parent is created). `constraints` and `objectives` are lists of specification mappings as in section 9; at least one of the two must be non-empty. `circular=True` treats the sequence as circular.
- **Output:** `status` (`solved`, `no_solution` or `no_solution_without_location`), `report_path` (`None` unless a report was written), `report_files`, `message`, `location` (`[start, end, strand]` or `None`, taken as DnaChisel reports it; in the example below its end, 11, lies beyond the 6-base sequence), `circular`, `operation` and `provenance`.
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.molecular_lab import write_dnachisel_no_solution_report

with tempfile.TemporaryDirectory() as folder:
    solvable = write_dnachisel_no_solution_report(
        "ATGGGTCTCATGCCCAAA",
        output_dir=Path(folder) / "solvable",
        constraints=[{"name": "AvoidPattern", "pattern": "GGTCTC"}],
        verbose=False,
    )
    print(solvable["status"], solvable["report_files"])
    stuck = write_dnachisel_no_solution_report(
        "GGTCTC",
        output_dir=Path(folder) / "stuck",
        constraints=[
            {"name": "AvoidPattern", "pattern": "GGTCTC"},
            {"name": "AvoidChanges", "location": [0, 6]},
        ],
        verbose=False,
    )
    print(stuck["status"], stuck["report_files"], stuck["location"])
```

- **Example output:**

```text
solved []
no_solution ['local_constraints_breaches.gb', 'logs.txt', 'plots.pdf'] [0, 11, 1]
```

## 20) `dnakit.molecular_lab.write_dnachisel_optimization_report()`

- **Purpose:** Optimize a sequence and write DnaChisel's optimization report folder.
- **Method:** Build the problem, resolve the constraints, optimize the objectives when there are any, and call DnaChisel's `write_optimization_report`. The folder holds the PDF summary, the before and after constraint and objective tables, the GenBank files and the logs. The two CSV tables are read back into the result.
- **API:** `write_dnachisel_optimization_report(sequence: 'str', *, output_dir: 'str | Path', constraints: 'Sequence[Mapping[str, Any]]' = (), objectives: 'Sequence[Mapping[str, Any]]' = (), circular: 'bool' = False, project_name: 'str' = 'DnaChisel optimization', verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the optional packages DnaChisel plus the DnaChisel report packages `pdf_reports` and `sequenticon` (extra `dnakit[bio-dnachisel-reports]`). They are not installed in the environment used to build this page, so the example shows the resulting `BackendUnavailableError`. Usage with the dependency installed is shown after the example.
- **Input:** `sequence` is an IUPAC DNA string. `output_dir` must not exist. At least one constraint or objective is required, given as specification mappings as in section 9. `circular` treats the sequence as circular. `project_name` is a non-empty string of at most 120 characters. If the constraints cannot be satisfied, `ConfigurationError` is raised; use `write_dnachisel_no_solution_report` for that case.
- **Output:** `report_path`, `report_files`, `sequence_length`, `edits`, `all_constraints_pass`, `constraints_before_after` and `objectives_before_after` (rows of the CSV tables), `circular`, `status`, `operation` and `provenance`.
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.exceptions import BackendUnavailableError
from dnakit.molecular_lab import write_dnachisel_optimization_report

with tempfile.TemporaryDirectory() as folder:
    try:
        write_dnachisel_optimization_report(
            "ATGGGTCTCATGCCCAAA",
            output_dir=Path(folder) / "report",
            constraints=[{"name": "AvoidPattern", "pattern": "GGTCTC"}],
            verbose=False,
        )
    except BackendUnavailableError as error:
        print(type(error).__name__, "pdf_reports could not be imported" in str(error))
```

- **Example output:**

```text
BackendUnavailableError True
```

Usage with the dependency installed:

```python
from dnakit.molecular_lab import write_dnachisel_optimization_report

result = write_dnachisel_optimization_report(
    "ATGGGTCTCATGCCCAAA",
    output_dir="dnachisel_report",
    constraints=[{"name": "AvoidPattern", "pattern": "GGTCTC"}],
    verbose=False,
)
print(result["edits"], result["all_constraints_pass"], result["report_files"])
```

## 21) `dnakit.molecular_lab.optimize_dnachisel_source_record()`

- **Purpose:** Optimize a record using the constraints and objectives annotated in it, and write the result.
- **Method:** Load the record with `DnaOptimizationProblem.from_record` (or `CircularDnaOptimizationProblem` with `circular=True`), which reads the constraints and objectives from the record's feature labels, as the DnaChisel command line does. A target ending in `.gb` receives an annotated GenBank record: the constraints are resolved, then the objectives optimized, then `to_record` writes the file. Any other target is a new folder or zip file written with `optimize_with_report`. The text summaries of the constraints and objectives are taken before and after.
- **API:** `optimize_dnachisel_source_record(source_path: 'str | Path', target_path: 'str | Path', *, circular: 'bool' = False, with_sequence_edits: 'bool' = False, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the optional package DnaChisel (extra `dnakit[bio-dnachisel]`). A non-`.gb` target (folder or zip report) also needs the report packages (extra `dnakit[bio-dnachisel-reports]`).
- **Input:** `source_path` is an existing GenBank (or other DnaChisel-readable) record whose features are labelled as DnaChisel specifications, for example `@AvoidPattern(GGTCTC)`. `target_path` must not exist and its parent folder must exist. `with_sequence_edits=True` annotates the edited positions in the GenBank output.
- **Output:** `source` (file name), `output_path`, `output_kind` (`genbank` or `report`), `circular`, `before` and `after` (each with `constraints` and `objectives` text summaries), `status`, `operation` and `provenance`.
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.molecular_lab import optimize_dnachisel_source_record

with tempfile.TemporaryDirectory() as folder:
    target = Path(folder) / "optimized.gb"
    result = optimize_dnachisel_source_record("annotated.gb", target, verbose=False)
    lines = target.read_text().splitlines()
print(result["output_kind"], result["source"])
print(result["before"]["constraints"].splitlines()[0])
print(result["after"]["constraints"].splitlines()[0])
print([line for line in lines if line.startswith("        1")][0])
```

- **Example output:**

```text
genbank annotated.gb
===> FAILURE: 1 constraints evaluations failed
===> SUCCESS - all constraints evaluations pass
        1 atgggtctaa tgcccaaagc tttaa
```

## 22) `dnakit.molecular_lab.dnacauldron_plot_leveled_graph()`

- **Purpose:** Draw a leveled assembly graph (nodes in columns by level, edges between them) into a PNG.
- **Method:** Call DnaCauldron's `AssemblyPlan.plot_leveled_graph` layout with a node-drawing callback that writes each node as a labelled rounded box, draw on an Agg matplotlib figure and save it as PNG (150 dpi, tight box). The node positions from the layout are returned as fractions of the plot area.
- **API:** `dnacauldron_plot_leveled_graph(output_path: 'str', levels: 'list[list[str]]', edges: 'list[list[str]]', *, width_factor: 'float' = 2.5, height_factor: 'float' = 2.0, scale: 'float' = 1.0, interlevel_shift: 'float' = 0.0, margin: 'float | None' = None, font_size: 'float' = 9.0, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the optional packages DnaCauldron (extra `dnakit[bio-dnacauldron]`) and matplotlib (extra `dnakit[bio-plots]`).
- **Input:** `output_path` must end in `.png` and not exist. `levels` is a list of 1 to 20 non-empty lists of unique node names (strings of 1 to 64 characters, at most 500 nodes in total). `edges` is a list of `[source, target]` pairs whose ends are nodes of `levels`. `width_factor` and `height_factor` are 0.1 to 50, `scale` 0.01 to 20, `interlevel_shift` -1 to 1, `margin` 0 to 1 or `None`, `font_size` 1 to 72.
- **Output:** `path` and `data` with `level_count`, `node_count`, `edge_count`, `levels`, `edges`, `positions` (node to `[x, y]`) and `provenance`.
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.molecular_lab import dnacauldron_plot_leveled_graph

levels = [["partA", "partB", "partC"], ["construct_1"]]
edges = [["partA", "construct_1"], ["partB", "construct_1"], ["partC", "construct_1"]]
with tempfile.TemporaryDirectory() as folder:
    result = dnacauldron_plot_leveled_graph(
        str(Path(folder) / "plan.png"), levels, edges, verbose=False
    )
    print(Path(result["path"]).read_bytes()[:4])
data = result["data"]
print(data["level_count"], data["node_count"], data["edge_count"])
print({node: [round(value, 3) for value in xy] for node, xy in data["positions"].items()})
```

- **Example output:**

```text
b'\x89PNG'
2 4 3
{'partA': [0.333, 0.25], 'partB': [0.333, 0.5], 'partC': [0.333, 0.75], 'construct_1': [0.667, 0.5]}
```

## 23) `dnakit.molecular_lab.dnacauldron_write_simulation_pdf_report()`

- **Purpose:** Write the DnaCauldron simulation summary as a PDF report.
- **Method:** Check the output path and the rows, import `pdf_reports`, build a pandas data frame with the columns `Outcome` and `Number of assemblies`, and call DnaCauldron's `reports.write_simulation_pdf_report`. These are the columns of the table that DnaCauldron's `AssemblyPlanSimulation` builds.
- **API:** `dnacauldron_write_simulation_pdf_report(output_path: 'str', simulation_info: 'list[dict[str, Any]]', *, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the optional package DnaCauldron (extra `dnakit[bio-dnacauldron]`) with its reports extra (`pip install "dnacauldron[reports]"`), which brings `pdf_reports` and, through it, WeasyPrint. These are not installed in the environment used to build this page, so the example shows the resulting `BackendUnavailableError`. Usage with the dependency installed is shown after the example.
- **Input:** `output_path` must end in `.pdf` and not exist. `simulation_info` is a list of 1 to 50 dictionaries with exactly the keys `Outcome` (a string of at most 64 characters, for example `Valid`, `Cancelled` or `Errored`) and `Number of assemblies` (a non-negative integer).
- **Output:** `path` and `data` with `rows`, `total_assemblies` and `provenance`.
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.exceptions import BackendUnavailableError
from dnakit.molecular_lab import dnacauldron_write_simulation_pdf_report

rows = [
    {"Outcome": "Valid", "Number of assemblies": 3},
    {"Outcome": "Errored", "Number of assemblies": 1},
]
with tempfile.TemporaryDirectory() as folder:
    try:
        dnacauldron_write_simulation_pdf_report(
            str(Path(folder) / "summary.pdf"), rows, verbose=False
        )
    except BackendUnavailableError as error:
        print(type(error).__name__, "pdf_reports could not be imported" in str(error))
```

- **Example output:**

```text
BackendUnavailableError True
```

Usage with the dependency installed:

```python
from dnakit.molecular_lab import dnacauldron_write_simulation_pdf_report

result = dnacauldron_write_simulation_pdf_report(
    "summary.pdf",
    [{"Outcome": "Valid", "Number of assemblies": 3}],
    verbose=False,
)
print(result["data"]["total_assemblies"])
```

## 24) `dnakit.molecular_lab.simulate_dnacauldron_plan()`

- **Purpose:** Simulate every assembly in a DnaCauldron plan and write the simulation report.
- **Method:** Import the part sequences from `sequences_path` into a `SequenceRepository`, import each sub-folder of `connectors_path` as a named collection, read the plan with `AssemblyPlan.from_spreadsheet`, run `simulate`, take `compute_stats()` and write the report with `write_report`, as the DnaCauldron command line does.
- **API:** `simulate_dnacauldron_plan(assembly_plan_path: 'str | Path', sequences_path: 'str | Path', target_path: 'str | Path', *, connectors_path: 'str | Path | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the optional package DnaCauldron (extra `dnakit[bio-dnacauldron]`). Writing the report folder worked here without `pdf_reports`.
- **Input:** `assembly_plan_path` is an existing CSV or spreadsheet file (at most 50 MB) with one row per assembly, in which a `class: ...` cell names the assembly class. `sequences_path` is an existing folder of part sequences. `connectors_path` is optional: a folder whose sub-folders each hold one connector collection. `target_path` must not exist, but its parent folder must; use a `.zip` name for a zip report.
- **Output:** `assembly_plan` (file name), `assemblies` (their number), `outcomes` (`valid_assemblies`, `cancelled_assemblies`, `errored_assemblies`), `report_path`, `status`, `operation` and `provenance`.
- **Example code:**

```python
import random
import tempfile
from pathlib import Path

from dnakit.molecular_lab import simulate_dnacauldron_plan

rng = random.Random(3)


def dna(length):
    return "".join(rng.choice("ACGT") for _ in range(length))


o0, o1, o2 = dna(40), dna(40), dna(40)  # shared overlaps between neighbouring parts
parts = {"partA": o2 + dna(160) + o0, "partB": o0 + dna(160) + o1, "partC": o1 + dna(160) + o2}
with tempfile.TemporaryDirectory() as folder:
    root = Path(folder)
    (root / "parts").mkdir()
    for name, sequence in parts.items():
        (root / "parts" / f"{name}.fa").write_text(f">{name}\n{sequence}\n")
    (root / "plan.csv").write_text("class: gibson_assembly,construct_1,partA,partB,partC\n")
    result = simulate_dnacauldron_plan(
        root / "plan.csv", root / "parts", root / "report", verbose=False
    )
    summary = root / "report" / "plan_simulation" / "plan_summary.csv"
    print(summary.is_file())
print(result["assembly_plan"], result["assemblies"])
print(result["outcomes"])
```

- **Example output:**

```text
True
plan.csv 1
{'cancelled_assemblies': 0, 'errored_assemblies': 0, 'valid_assemblies': 1}
```
