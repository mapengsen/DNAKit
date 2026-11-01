# DendroPy and khmer: tree formats, PAUP*, and k-mer tables

This group of functions is a set of adapters over two optional packages. [DendroPy](https://github.com/jeetsukumaran/DendroPy) is a Python library for phylogenetic computing; `dnakit.phylogeny` uses it to convert sequence and tree files between FASTA, Newick, NEXUS, NeXML and PHYLIP, to run its `dendropy-format`, `sumtrees` and `sumlabels` console applications, to fetch records from NCBI Entrez, and to drive the external programs MUSCLE (alignment) and PAUP* (tree and model estimation, pruning, parsimony scores). [khmer](https://github.com/dib-lab/khmer) is a library for approximate k-mer counting; `dnakit.microbiome` uses it to size khmer tables and to build and save counting (countgraph) and presence (nodegraph) tables from reads. The functions are imported with `from dnakit.phylogeny import ...` and `from dnakit.microbiome import ...`.

DendroPy and khmer are BSD-licensed according to their package metadata. MUSCLE and PAUP* are separate programs that DNAKit does not distribute; their licences are not checked here.

The dependencies are optional and are imported only when the matching function is called. If a package is missing, the function raises an error that includes the install command (`dnakit[bio-dendropy]` or `dnakit[bio-khmer]`). Output files are always new: an existing output path is refused, never overwritten.

This environment has DendroPy 5.1.0 (with the `dendropy-format`, `sumtrees` and `sumlabels` scripts) and khmer 2.1.1, but neither MUSCLE nor PAUP*, and no network access. Twelve examples therefore run the real computation (sections 2, 4 to 8 and 14 to 19). Sections 1, 3 and 9 to 13 show the real error that the function raises when the program or the input is missing, and each of them is followed by an unexecuted "Usage with the dependency installed" block.

## Example data

The examples use five tiny declared files: a four-sequence DNA alignment, a two-sequence file with ambiguity codes, a Newick file with three trees, a Newick target tree and a FASTA file of five short reads (two identical, one short, one with an `N`). Every output is written to a temporary directory.

**Example file:** `aln.fa`

```text
>t1
ACGTACGTAC
>t2
ACGTACGAAC
>t3
ACCTACGTTC
>t4
TCCTACGTTC
```

**Example file:** `iupac.fa`

```text
>a
ACGTRYAC
>b
ACGTACNC
```

**Example file:** `trees.nwk`

```text
((t1:0.1,t2:0.1)0.9:0.2,(t3:0.1,t4:0.1)0.8:0.2);
((t1:0.1,t2:0.1)0.7:0.2,(t3:0.1,t4:0.1)0.6:0.2);
((t1:0.1,t3:0.1)0.5:0.2,(t2:0.1,t4:0.1)0.4:0.2);
```

**Example file:** `target.nwk`

```text
((t1,t2),(t3,t4));
```

**Example file:** `reads.fa`

```text
>r1
ACGTACGTACGTAGCT
>r2
ACGTACGTACGTAGCT
>r3
TTGCAGGTCCAT
>r4_short
ACG
>r5_n
ACGTNCGTACGT
```

## 1) `dnakit.phylogeny.align_with_muscle()`

- **Purpose:** Align DNA sequences with the external MUSCLE program and save the alignment as FASTA.
- **Method:** Reads the FASTA file into a dendropy `DnaCharacterMatrix`, starts MUSCLE through dendropy's `muscle_align` (the sequences go in on standard input and the alignment is read back from MUSCLE's standard output), and writes the aligned matrix as FASTA to a new file. The executable is looked up from `muscle_path`, then on PATH, then in the folder of the running interpreter. A non-zero exit raises `BackendExecutionError` carrying MUSCLE's message.
- **API:** `align_with_muscle(fasta: 'str', output: 'str', *, muscle_args: 'Sequence[str] | None' = None, muscle_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** needs the optional package dendropy 5.1 (extra `dnakit[bio-dendropy]`) and the external program `muscle` on PATH (or passed as `muscle_path`). MUSCLE is not installed in the environment used to build this page, so the example below shows the real error that is raised when the program cannot be found; a call with the program installed follows it (not executed).
- **Input:** `fasta` is an existing FASTA file of DNA sequences; `output` is a path that does not exist yet in an existing folder; `muscle_args` is an optional list of non-empty strings added to the MUSCLE command line (none by default, which is the MUSCLE 3.8 style); `muscle_path` is an optional path to the executable.
- **Output:** A dictionary with `status`, `backend`, `version`, `command`, `sequence_count`, `alignment_length` and `outputs` (the new FASTA file).
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.exceptions import BackendUnavailableError
from dnakit.phylogeny import align_with_muscle

with tempfile.TemporaryDirectory() as folder:
    try:
        align_with_muscle("aln.fa", str(Path(folder) / "aligned.fa"), muscle_path="no_such_dir/muscle", verbose=False)
    except BackendUnavailableError as error:
        print(type(error).__name__, str(error).split(":")[0])
```

- **Example output:**

```text
BackendUnavailableError [BACKEND_UNAVAILABLE] muscle executable does not exist
```

Usage with the dependency installed:

```python
from dnakit.phylogeny import align_with_muscle

result = align_with_muscle("aln.fa", "aligned.fa", verbose=False)
print(result["sequence_count"], result["alignment_length"])
```

## 2) `dnakit.phylogeny.convert_format()`

- **Purpose:** Convert a DNA sequence file or a tree file from one format to another.
- **Method:** Follows the steps of the `dendropy-format` application but stays in the current process: reads the file with `dendropy.DataSet.get`, optionally recodes ambiguous and polymorphic DNA codes, and writes it with `DataSet.write`. It deliberately does not call the application's `convert()` function, which closes `sys.stdout`. With `recode_uncertain="missing"` the codes become `?` and with `"gap"` they become `-`.
- **API:** `convert_format(source: 'str', output: 'str', *, input_format: 'str', output_format: 'str | None' = None, data_type: 'str | None' = None, unquoted_underscores: 'bool' = False, recode_uncertain: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** needs the optional package dendropy 5.1 (extra `dnakit[bio-dendropy]`). Works after plain pip install of that extra; no external program is needed.
- **Input:** `source` is an existing file; `output` is a new path; `input_format` is one of `fasta`, `newick`, `nexml`, `nexus` or a `phylip-*` name (`phylip-strict`, `phylip-relaxed-singlespace`, `phylip-relaxed-multispace` and the three `-interleaved` variants); `output_format` is one of `fasta`, `newick`, `nexml`, `nexus`, `phylip`, `phylip-strict` and defaults to the reader's format; `data_type` (`dna`, `rna` or `standard`) is required for FASTA and PHYLIP input; `unquoted_underscores` applies to NEXUS and Newick output; `recode_uncertain` is `missing` or `gap`.
- **Output:** A dictionary with `status`, `backend`, `version`, `operation`, `input_format`, `output_format`, `data_type`, `tree_lists` and `character_matrices` (the numbers read) and `outputs`.
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.phylogeny import convert_format

with tempfile.TemporaryDirectory() as folder:
    output = Path(folder) / "iupac.nex"
    result = convert_format(
        "iupac.fa", str(output), input_format="fasta", output_format="nexus", data_type="dna",
        recode_uncertain="missing", verbose=False,
    )
    print(result["status"], result["input_format"], result["output_format"], result["character_matrices"])
    for line in output.read_text().splitlines():
        if len(line.split()) == 2 and line.split()[0] in ("a", "b"):
            print(" ".join(line.split()))
```

- **Example output:**

```text
completed fasta nexus 1
a ACGT??AC
b ACGTAC?C
```

## 3) `dnakit.phylogeny.fetch_entrez_records()`

- **Purpose:** Fetch sequence records from an NCBI database by numeric Entrez UID.
- **Method:** Validates the arguments, then calls `dendropy.interop.entrez.efetch`, which sends an E-utilities efetch request to NCBI over HTTP, reads the whole response and closes it. The identifiers are joined with commas. Nothing is cached and nothing is written to disk. `record_count` counts the `>` header lines and is only set when `rettype` is `fasta` and `retmode` is `text`.
- **API:** `fetch_entrez_records(database: 'str', ids: 'Sequence[str | int]', *, rettype: 'str' = 'fasta', retmode: 'str' = 'text', email: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** needs internet access to NCBI (no key is passed; `email` is optional) and the optional package dendropy 5.1 (extra `dnakit[bio-dendropy]`). The build environment is offline, so the example shows the real input check that happens before any request; a call with network access follows it (not executed).
- **Input:** `database` is a non-empty name that does not start with `-` (for example `nuccore`); `ids` is a list of 1 to 200 numeric UIDs (integers or digit strings; accessions such as `NM_000546` are refused); `rettype` is a non-empty name (default `fasta`); `retmode` is `text` or `xml`; `email` is optional and must contain `@`.
- **Output:** A dictionary with `status`, `backend`, `service`, `database`, `rettype`, `retmode`, `identifier_count`, `record_count`, `access_date` (UTC date), `url` (the final request URL) and `records` (the response text).
- **Example code:**

```python
from dnakit.exceptions import ConfigurationError
from dnakit.phylogeny import fetch_entrez_records

try:
    fetch_entrez_records("nuccore", ["NM_000546"], verbose=False)
except ConfigurationError as error:
    print(type(error).__name__, str(error).split(". ")[0])
```

- **Example output:**

```text
ConfigurationError [CONFIGURATION_ERROR] Entrez identifiers are numeric UIDs; got 'NM_000546'
```

Usage with the dependency installed:

```python
from dnakit.phylogeny import fetch_entrez_records

# 1234567 stands for any numeric nuccore UID, for example one found with an ESearch query.
result = fetch_entrez_records("nuccore", [1234567], email="you@example.org", verbose=False)
print(result["record_count"], result["access_date"])
print(result["records"][:80])
```

## 4) `dnakit.phylogeny.format_to_fasta()`

- **Purpose:** Convert a sequence file (or the sequences of a NEXUS file) to FASTA.
- **Method:** Calls `convert_format` with `output_format="fasta"`; this is `dendropy-format --to fasta`.
- **API:** `format_to_fasta(source: 'str', output: 'str', *, input_format: 'str', data_type: 'str | None' = None, recode_uncertain: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** needs the optional package dendropy 5.1 (extra `dnakit[bio-dendropy]`). Works after installing that extra; no external program is needed.
- **Input:** `source`, `output`, `input_format`, `data_type` and `recode_uncertain` as in `convert_format` (FASTA and PHYLIP input need `data_type`).
- **Output:** The same dictionary as `convert_format`, with `output_format` equal to `fasta`.
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.phylogeny import format_to_fasta

with tempfile.TemporaryDirectory() as folder:
    output = Path(folder) / "gapped.fa"
    result = format_to_fasta(
        "iupac.fa", str(output), input_format="fasta", data_type="dna", recode_uncertain="gap", verbose=False
    )
    print(result["output_format"])
    print(output.read_text().strip())
```

- **Example output:**

```text
fasta
>a
ACGT--AC
>b
ACGTAC-C
```

## 5) `dnakit.phylogeny.format_to_newick()`

- **Purpose:** Convert a tree file (or the trees of a NEXUS file) to Newick.
- **Method:** Calls `convert_format` with `output_format="newick"`; this is `dendropy-format --to newick`. Every tree is written on its own line.
- **API:** `format_to_newick(source: 'str', output: 'str', *, input_format: 'str', data_type: 'str | None' = None, unquoted_underscores: 'bool' = False, recode_uncertain: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** needs the optional package dendropy 5.1 (extra `dnakit[bio-dendropy]`). Works after installing that extra; no external program is needed.
- **Input:** `source`, `output`, `input_format`, `data_type`, `unquoted_underscores` and `recode_uncertain` as in `convert_format`.
- **Output:** The same dictionary as `convert_format`, with `output_format` equal to `newick`.
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.phylogeny import format_to_newick

with tempfile.TemporaryDirectory() as folder:
    output = Path(folder) / "trees_out.nwk"
    result = format_to_newick("trees.nwk", str(output), input_format="newick", verbose=False)
    print(result["output_format"], result["tree_lists"], result["character_matrices"])
    print(output.read_text().splitlines()[0])
```

- **Example output:**

```text
newick 1 0
((t1:0.1,t2:0.1)0.9:0.2,(t3:0.1,t4:0.1)0.8:0.2);
```

## 6) `dnakit.phylogeny.format_to_nexml()`

- **Purpose:** Convert a sequence or tree file to NeXML.
- **Method:** Calls `convert_format` with `output_format="nexml"`; this is `dendropy-format --to nexml`. There is no `unquoted_underscores` option because it applies only to NEXUS and Newick.
- **API:** `format_to_nexml(source: 'str', output: 'str', *, input_format: 'str', data_type: 'str | None' = None, recode_uncertain: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** needs the optional package dendropy 5.1 (extra `dnakit[bio-dendropy]`). Works after installing that extra; no external program is needed.
- **Input:** `source`, `output`, `input_format`, `data_type` and `recode_uncertain` as in `convert_format`.
- **Output:** The same dictionary as `convert_format`, with `output_format` equal to `nexml`.
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.phylogeny import format_to_nexml

with tempfile.TemporaryDirectory() as folder:
    output = Path(folder) / "aln.xml"
    result = format_to_nexml("aln.fa", str(output), input_format="fasta", data_type="dna", verbose=False)
    text = output.read_text()
    print(result["output_format"], text.startswith("<?xml"), text.count("<otu "), text.count("<row"))
```

- **Example output:**

```text
nexml True 4 4
```

## 7) `dnakit.phylogeny.format_to_nexus()`

- **Purpose:** Convert a sequence or tree file to NEXUS.
- **Method:** Calls `convert_format` with `output_format="nexus"`; this is `dendropy-format --to nexus`. The output has a TAXA block and a TREES or CHARACTERS block.
- **API:** `format_to_nexus(source: 'str', output: 'str', *, input_format: 'str', data_type: 'str | None' = None, unquoted_underscores: 'bool' = False, recode_uncertain: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** needs the optional package dendropy 5.1 (extra `dnakit[bio-dendropy]`). Works after installing that extra; no external program is needed.
- **Input:** `source`, `output`, `input_format`, `data_type`, `unquoted_underscores` and `recode_uncertain` as in `convert_format`.
- **Output:** The same dictionary as `convert_format`, with `output_format` equal to `nexus`.
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.phylogeny import format_to_nexus

with tempfile.TemporaryDirectory() as folder:
    output = Path(folder) / "trees.nex"
    result = format_to_nexus("trees.nwk", str(output), input_format="newick", verbose=False)
    print(result["output_format"], result["tree_lists"])
    for line in output.read_text().splitlines():
        if line.strip().startswith("TREE"):
            print(line.strip())
```

- **Example output:**

```text
nexus 1
TREE 1 = ((t1:0.1,t2:0.1)0.9:0.2,(t3:0.1,t4:0.1)0.8:0.2);
TREE 2 = ((t1:0.1,t2:0.1)0.7:0.2,(t3:0.1,t4:0.1)0.6:0.2);
TREE 3 = ((t1:0.1,t3:0.1)0.5:0.2,(t2:0.1,t4:0.1)0.4:0.2);
```

## 8) `dnakit.phylogeny.format_to_phylip()`

- **Purpose:** Convert a sequence file to relaxed PHYLIP.
- **Method:** Calls `convert_format` with `output_format="phylip"` (dendropy's relaxed PHYLIP writer); this is `dendropy-format --to phylip`. The strict variant is available in `convert_format` as `output_format="phylip-strict"`.
- **API:** `format_to_phylip(source: 'str', output: 'str', *, input_format: 'str', data_type: 'str | None' = None, recode_uncertain: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** needs the optional package dendropy 5.1 (extra `dnakit[bio-dendropy]`). Works after installing that extra; no external program is needed.
- **Input:** `source`, `output`, `input_format`, `data_type` and `recode_uncertain` as in `convert_format`.
- **Output:** The same dictionary as `convert_format`, with `output_format` equal to `phylip`.
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.phylogeny import format_to_phylip

with tempfile.TemporaryDirectory() as folder:
    output = Path(folder) / "aln.phy"
    result = format_to_phylip("aln.fa", str(output), input_format="fasta", data_type="dna", verbose=False)
    print(result["output_format"])
    print(output.read_text().strip())
```

- **Example output:**

```text
phylip
4 10
t1  ACGTACGTAC
t2  ACGTACGAAC
t3  ACCTACGTTC
t4  TCCTACGTTC
```

## 9) `dnakit.phylogeny.paup_estimate_model()`

- **Purpose:** Estimate the substitution-model parameters of a DNA alignment with PAUP*.
- **Method:** Reads the alignment (FASTA) and, if `tree_file` is given, its first tree (the file is taken as NEXUS when its first non-empty line starts with `#NEXUS`, otherwise as Newick). Calls dendropy's `paup.estimate_model`, which fits the model on that tree, or on a tree first estimated with `criterion`, and parses PAUP*'s report. The tree PAUP* used is written as Newick.
- **API:** `paup_estimate_model(alignment: 'str', output_tree: 'str', *, tree_file: 'str | None' = None, criterion: 'str' = 'likelihood', substitution_types: 'int' = 6, unequal_base_freqs: 'bool' = True, gamma_rates: 'bool' = True, prop_invar: 'bool' = True, user_branch_lengths: 'bool' = True, paup_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** needs the optional package dendropy 5.1 (extra `dnakit[bio-dendropy]`) and the external program `paup` (PAUP* 4.0) on PATH, or passed as `paup_path`. PAUP* is not installed in the environment used to build this page, so the example shows the real error raised when the program cannot be found; a call with PAUP* installed follows it (not executed). The executable path may contain only letters, digits and `._/+:@%-`, because dendropy starts PAUP* through a shell in some of its functions.
- **Input:** `alignment` is an existing DNA FASTA file; `output_tree` is a new path; `tree_file` is an optional Newick or NEXUS file; `criterion` is `likelihood`, `parsimony`, `nj` or `upgma`; `substitution_types` is 1, 2 or 6; `unequal_base_freqs`, `gamma_rates`, `prop_invar` and `user_branch_lengths` are booleans; `paup_path` is optional.
- **Output:** A dictionary with `status`, `backend`, `operation`, `version`, `outputs` (the Newick file), `substitution_types` and `estimates` (`log_likelihood`, `kappa`, `rAC` to `rGT`, `pA` to `pT`, `prop_invar`, `alpha`; each a float or `None` when PAUP* did not report it).
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.exceptions import BackendUnavailableError
from dnakit.phylogeny import paup_estimate_model

with tempfile.TemporaryDirectory() as folder:
    try:
        paup_estimate_model("aln.fa", str(Path(folder) / "model_tree.nwk"), paup_path="no_such_dir/paup", verbose=False)
    except BackendUnavailableError as error:
        print(type(error).__name__, str(error).split(":")[0])
```

- **Example output:**

```text
BackendUnavailableError [BACKEND_UNAVAILABLE] paup executable does not exist
```

Usage with the dependency installed:

```python
from dnakit.phylogeny import paup_estimate_model

result = paup_estimate_model("aln.fa", "model_tree.nwk", tree_file="start_tree.nwk", substitution_types=2, verbose=False)
print(result["estimates"]["log_likelihood"], result["estimates"]["kappa"])
```

## 10) `dnakit.phylogeny.paup_estimate_tree()`

- **Purpose:** Estimate a tree from a DNA alignment with PAUP* and save it as Newick.
- **Method:** Reads the alignment (FASTA) and calls dendropy's `paup.estimate_tree` with the criterion and the model options, then writes the tree as Newick. The model options (`substitution_types`, `unequal_base_freqs`, `gamma_rates`, `prop_invar`) apply to likelihood estimation. The program is looked up before the alignment is read, and a missing PAUP* raises `BackendUnavailableError`.
- **API:** `paup_estimate_tree(alignment: 'str', output_tree: 'str', *, criterion: 'str' = 'likelihood', substitution_types: 'int' = 6, unequal_base_freqs: 'bool' = True, gamma_rates: 'bool' = True, prop_invar: 'bool' = True, paup_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** needs the optional package dendropy 5.1 (extra `dnakit[bio-dendropy]`) and the external program `paup` (PAUP* 4.0) on PATH, or passed as `paup_path`. PAUP* is not installed in the environment used to build this page, so the example shows the real error raised when the program cannot be found; a call with PAUP* installed follows it (not executed). The executable path may contain only letters, digits and `._/+:@%-`, because dendropy starts PAUP* through a shell in some of its functions.
- **Input:** `alignment` is an existing DNA FASTA file; `output_tree` is a new path; `criterion` is `likelihood` (default), `parsimony`, `nj` or `upgma`; `substitution_types` is 1, 2 or 6 (PAUP*'s `lset nst`); the three model switches are booleans; `paup_path` is optional.
- **Output:** A dictionary with `status`, `backend`, `operation`, `version`, `outputs` (the Newick file), `criterion`, `substitution_types` and `taxon_count`.
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.exceptions import BackendUnavailableError
from dnakit.phylogeny import paup_estimate_tree

with tempfile.TemporaryDirectory() as folder:
    try:
        paup_estimate_tree("aln.fa", str(Path(folder) / "ml_tree.nwk"), paup_path="no_such_dir/paup", verbose=False)
    except BackendUnavailableError as error:
        print(type(error).__name__, str(error).split(":")[0])
```

- **Example output:**

```text
BackendUnavailableError [BACKEND_UNAVAILABLE] paup executable does not exist
```

Usage with the dependency installed:

```python
from dnakit.phylogeny import paup_estimate_tree

result = paup_estimate_tree("aln.fa", "nj_tree.nwk", criterion="nj", verbose=False)
print(result["criterion"], result["taxon_count"], result["outputs"])
```

## 11) `dnakit.phylogeny.paup_estimate_ultrametric_tree()`

- **Purpose:** Estimate an ultrametric (clock-constrained) tree with PAUP*.
- **Method:** Reads the alignment (FASTA) and, if `topology_tree` is given, its first tree. Calls dendropy's `paup.estimate_ultrametric_tree`: without a topology it builds a neighbour-joining tree first and then estimates clock-constrained branch lengths; with a topology it keeps that topology and estimates only the branch lengths. The tree is written as Newick with branch lengths.
- **API:** `paup_estimate_ultrametric_tree(alignment: 'str', output_tree: 'str', *, topology_tree: 'str | None' = None, paup_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** needs the optional package dendropy 5.1 (extra `dnakit[bio-dendropy]`) and the external program `paup` (PAUP* 4.0) on PATH, or passed as `paup_path`. PAUP* is not installed in the environment used to build this page, so the example shows the real error raised when the program cannot be found; a call with PAUP* installed follows it (not executed). The executable path may contain only letters, digits and `._/+:@%-`, because dendropy starts PAUP* through a shell in some of its functions.
- **Input:** `alignment` is an existing DNA FASTA file; `output_tree` is a new path; `topology_tree` is an optional Newick or NEXUS file whose first tree is used; `paup_path` is optional.
- **Output:** A dictionary with `status`, `backend`, `operation`, `version`, `outputs` (the Newick file), `taxon_count` and `topology_given`.
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.exceptions import BackendUnavailableError
from dnakit.phylogeny import paup_estimate_ultrametric_tree

with tempfile.TemporaryDirectory() as folder:
    try:
        paup_estimate_ultrametric_tree(
            "aln.fa", str(Path(folder) / "clock_tree.nwk"), paup_path="no_such_dir/paup", verbose=False
        )
    except BackendUnavailableError as error:
        print(type(error).__name__, str(error).split(":")[0])
```

- **Example output:**

```text
BackendUnavailableError [BACKEND_UNAVAILABLE] paup executable does not exist
```

Usage with the dependency installed:

```python
from dnakit.phylogeny import paup_estimate_ultrametric_tree

result = paup_estimate_ultrametric_tree("aln.fa", "clock_tree.nwk", topology_tree="target.nwk", verbose=False)
print(result["topology_given"], result["taxon_count"])
```

## 12) `dnakit.phylogeny.paup_prune_taxa()`

- **Purpose:** Delete named taxa from every tree of a file with PAUP* and save the pruned trees as NEXUS.
- **Method:** Reads all trees of the file into one taxon namespace and looks up every name in `taxa` (an unknown name raises `ConfigurationError` before PAUP* is looked up). Then calls dendropy's `paup.prune_taxa_from_trees` and writes the pruned trees as NEXUS with a taxa block.
- **API:** `paup_prune_taxa(trees_file: 'str', taxa: 'Sequence[str]', output_trees: 'str', *, paup_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** needs the optional package dendropy 5.1 (extra `dnakit[bio-dendropy]`) and the external program `paup` (PAUP* 4.0) on PATH, or passed as `paup_path`. PAUP* is not installed in the environment used to build this page, so the example shows the real error raised when the program cannot be found; a call with PAUP* installed follows it (not executed). The executable path may contain only letters, digits and `._/+:@%-`, because dendropy starts PAUP* through a shell in some of its functions.
- **Input:** `trees_file` is an existing Newick or NEXUS file; `taxa` is a non-empty list of single-line labels that all occur in the trees; `output_trees` is a new path; `paup_path` is optional.
- **Output:** A dictionary with `status`, `backend`, `operation`, `version`, `outputs` (the NEXUS file), `removed` (the labels) and `tree_count`.
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.exceptions import BackendUnavailableError, ConfigurationError
from dnakit.phylogeny import paup_prune_taxa

with tempfile.TemporaryDirectory() as folder:
    output = str(Path(folder) / "pruned.nex")
    try:
        paup_prune_taxa("trees.nwk", ["t4"], output, paup_path="no_such_dir/paup", verbose=False)
    except BackendUnavailableError as error:
        print(type(error).__name__, str(error).split(":")[0])
    try:
        paup_prune_taxa("trees.nwk", ["no_such_taxon"], output, paup_path="no_such_dir/paup", verbose=False)
    except ConfigurationError as error:
        print(type(error).__name__, str(error).split(":")[0])
```

- **Example output:**

```text
BackendUnavailableError [BACKEND_UNAVAILABLE] paup executable does not exist
ConfigurationError [CONFIGURATION_ERROR] Taxon not found in the trees
```

Usage with the dependency installed:

```python
from dnakit.phylogeny import paup_prune_taxa

result = paup_prune_taxa("trees.nwk", ["t4"], "pruned.nex", verbose=False)
print(result["removed"], result["tree_count"])
```

## 13) `dnakit.phylogeny.paup_pscore_trees()`

- **Purpose:** Compute the parsimony score (PAUP* `pscore`) of every tree in a file against a DNA alignment.
- **Method:** Reads the trees and the alignment into one taxon namespace and calls dendropy's `paup.pscore_trees`; the scores are converted to integers. Each `pset_options` item is written as part of a `pset` command, so it must match the plain token pattern `[A-Za-z0-9_.=:+-]+` (for example `gapmode=missing`); anything else is refused before PAUP* is looked up. dendropy's own `pscore_option_list` is not exposed because it joins its options to the score-file path without a separator.
- **API:** `paup_pscore_trees(trees_file: 'str', alignment: 'str', *, pset_options: 'Sequence[str] | None' = None, paup_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** needs the optional package dendropy 5.1 (extra `dnakit[bio-dendropy]`) and the external program `paup` (PAUP* 4.0) on PATH, or passed as `paup_path`. PAUP* is not installed in the environment used to build this page, so the example shows the real error raised when the program cannot be found; a call with PAUP* installed follows it (not executed). The executable path may contain only letters, digits and `._/+:@%-`, because dendropy starts PAUP* through a shell in some of its functions. Invalid `pset_options` are refused without PAUP*, which the example also shows.
- **Input:** `trees_file` is an existing Newick or NEXUS file; `alignment` is an existing DNA FASTA file; `pset_options` is an optional list of plain tokens; `paup_path` is optional.
- **Output:** A dictionary with `status`, `backend`, `operation`, `version`, `outputs` (always empty), `scores` (one integer per tree) and `tree_count`.
- **Example code:**

```python
from dnakit.exceptions import BackendUnavailableError, ConfigurationError
from dnakit.phylogeny import paup_pscore_trees

try:
    paup_pscore_trees("trees.nwk", "aln.fa", pset_options=["gapmode=missing"], paup_path="no_such_dir/paup", verbose=False)
except BackendUnavailableError as error:
    print(type(error).__name__, str(error).split(":")[0])
try:
    paup_pscore_trees("trees.nwk", "aln.fa", pset_options=["gapmode=missing; quit"], paup_path="no_such_dir/paup", verbose=False)
except ConfigurationError as error:
    print(type(error).__name__, str(error).split(";")[0])
```

- **Example output:**

```text
BackendUnavailableError [BACKEND_UNAVAILABLE] paup executable does not exist
ConfigurationError [CONFIGURATION_ERROR] pset_options must be plain PAUP* tokens such as gapmode=missing
```

Usage with the dependency installed:

```python
from dnakit.phylogeny import paup_pscore_trees

result = paup_pscore_trees("trees.nwk", "aln.fa", pset_options=["gapmode=missing"], verbose=False)
print(result["tree_count"], result["scores"])
```

## 14) `dnakit.phylogeny.run_dendropy_format()`

- **Purpose:** Run the `dendropy-format` console application on a file and save its standard output as a new file.
- **Method:** Builds a fixed argument list (`--from`, `--to`, `--data-type`, `--unquoted-underscores`, `--recode-uncertain`, then the source path) and runs `dendropy-format` without a shell through DNAKit's bounded command runner. If the exit code is 0 the captured output is written to the new file; otherwise no file is written and `status` is `failed`. `recode_uncertain="missing"` is refused with `ConfigurationError` because `dendropy-format` 5.1 crashes with it; use `"gap"`, or `convert_format`, which writes `?`. The same conversion without a subprocess is `convert_format`.
- **API:** `run_dendropy_format(source: 'str', output: 'str', *, input_format: 'str', output_format: 'str | None' = None, data_type: 'str | None' = None, unquoted_underscores: 'bool' = False, recode_uncertain: 'str | None' = None, dendropy_format_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** needs the optional package dendropy 5.1 (extra `dnakit[bio-dendropy]`), which installs the `dendropy-format` console script; the script is found on PATH or in the folder of the running interpreter, or from `dendropy_format_path`.
- **Input:** `source` is an existing file; `output` is a new path; `input_format`, `output_format`, `data_type` and `unquoted_underscores` as in `convert_format` (`data_type` is required for FASTA and PHYLIP input); `recode_uncertain` is `gap`; `dendropy_format_path` is optional.
- **Output:** A dictionary with `status` (`completed` or `failed`), `backend`, `version`, `command`, `return_code`, `outputs` (the new file, or an empty list on failure) and `log_tail` (empty on success, the end of the program output on failure).
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.phylogeny import run_dendropy_format

with tempfile.TemporaryDirectory() as folder:
    output = Path(folder) / "trees.nex"
    result = run_dendropy_format("trees.nwk", str(output), input_format="newick", output_format="nexus", verbose=False)
    print(result["status"], result["backend"], result["return_code"], len(result["outputs"]))
    for line in output.read_text().splitlines():
        if line.strip().startswith("TREE"):
            print(line.strip())
```

- **Example output:**

```text
completed dendropy-format 0 1
TREE 1 = ((t1:0.1,t2:0.1)0.9:0.2,(t3:0.1,t4:0.1)0.8:0.2);
TREE 2 = ((t1:0.1,t2:0.1)0.7:0.2,(t3:0.1,t4:0.1)0.6:0.2);
TREE 3 = ((t1:0.1,t3:0.1)0.5:0.2,(t2:0.1,t4:0.1)0.4:0.2);
```

## 15) `dnakit.phylogeny.summarize_tree_samples()`

- **Purpose:** Summarise a sample of trees (for example MCMC output) into one summary tree with support values.
- **Method:** Builds a fixed argument list for the `sumtrees` application (`--quiet`, `--burnin`, `--summary-target` or `--target-tree-filepath`, `--labels`, `--decimals`, `--output`, then the input files) and runs it without a shell. `summary_target` is `consensus` (default), `mcct` or `msct`; `target_tree` instead scores the clades of one given tree, and giving both is refused with `ConfigurationError` because sumtrees crashes on that combination. The summary is written by sumtrees to the new output file.
- **API:** `summarize_tree_samples(tree_files: 'Sequence[str]', output: 'str', *, input_format: 'str | None' = None, burnin: 'int' = 0, rooted: 'bool | None' = None, weighted_trees: 'bool' = False, summary_target: 'str | None' = None, target_tree: 'str | None' = None, min_clade_freq: 'float | None' = None, edges: 'str | None' = None, labels: 'str' = 'support', percentages: 'bool' = False, decimals: 'int' = 8, output_format: 'str | None' = None, sumtrees_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** needs the optional package dendropy 5.1 (extra `dnakit[bio-dendropy]`), which installs the `sumtrees` console script; it is found on PATH or in the folder of the running interpreter, or from `sumtrees_path`.
- **Input:** `tree_files` is a non-empty list of existing tree files; `output` is a new path; `input_format` is `nexus/newick`, `nexus`, `newick`, `phylip` or `nexml`; `burnin` is a non-negative integer (trees dropped from the start of each file); `rooted` is `True` (force rooted), `False` (force unrooted) or `None`; `weighted_trees`, `percentages` are booleans; `summary_target` or `target_tree`; `min_clade_freq` is between 0 and 1; `edges` is `mean-length`, `median-length`, `mean-age`, `median-age`, `support`, `keep` or `clear`; `labels` is `support`, `keep` or `clear`; `decimals` is a non-negative integer; `output_format` is `nexus`, `newick`, `phylip` or `nexml`; `sumtrees_path` is optional.
- **Output:** A dictionary with `status` (`completed` when the exit code is 0 and the output file exists, otherwise `failed`), `backend`, `version`, `command`, `return_code`, `outputs` and `log_tail`.
- **Example code:**

```python
import re
import tempfile
from pathlib import Path

from dnakit.phylogeny import summarize_tree_samples

with tempfile.TemporaryDirectory() as folder:
    output = Path(folder) / "summary.nwk"
    result = summarize_tree_samples(
        ["trees.nwk"], str(output), input_format="newick", rooted=True, output_format="newick", verbose=False
    )
    print(result["status"], result["backend"], result["return_code"])
    print(re.findall(r"\)([0-9.]+):", output.read_text()))
```

- **Example output:**

```text
completed sumtrees 0
['0.66666667', '0.66666667', '1.00000000']
```

## 16) `dnakit.phylogeny.transfer_tree_labels()`

- **Purpose:** Copy the labels that source trees carry on matching splits onto one target tree.
- **Method:** Builds a fixed argument list for the `sumlabels` application (`--quiet`, `--target`, `--separator=...`, optional flags, `--output`, then the source files) and runs it without a shell. Each split of the target tree receives the labels found on the same split in the source trees, joined by `separator`. The output is NEXUS with a taxa block by default, or Newick with `to_newick=True`.
- **API:** `transfer_tree_labels(tree_files: 'Sequence[str]', target_tree: 'str', output: 'str', *, separator: 'str' = '/', preserve_target_labels: 'bool' = False, rooted: 'bool | None' = None, include_taxa_block: 'bool' = True, additional_comments: 'str | None' = None, to_newick: 'bool' = False, sumlabels_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** needs the optional package dendropy 5.1 (extra `dnakit[bio-dendropy]`), which installs the `sumlabels` console script; it is found on PATH or in the folder of the running interpreter, or from `sumlabels_path`.
- **Input:** `tree_files` is a non-empty list of existing tree files with labelled internal nodes; `target_tree` is an existing tree file; `output` is a new path; `separator` is a non-empty single-line string (default `/`); `preserve_target_labels`, `include_taxa_block` and `to_newick` are booleans; `rooted` is `True`, `False` or `None`; `additional_comments` is an optional single-line string; `sumlabels_path` is optional.
- **Output:** A dictionary with `status` (`completed` when the exit code is 0 and the output file exists, otherwise `failed`), `backend`, `version`, `command`, `return_code`, `outputs` and `log_tail`. `version` can be `None` because `sumlabels --version` does not print a plain version number.
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.phylogeny import transfer_tree_labels

with tempfile.TemporaryDirectory() as folder:
    output = Path(folder) / "labelled.nwk"
    result = transfer_tree_labels(["trees.nwk"], "target.nwk", str(output), rooted=True, to_newick=True, verbose=False)
    print(result["status"], result["backend"], result["return_code"])
    print(output.read_text().strip())
```

- **Example output:**

```text
completed sumlabels 0
[&R] ((t1,t2)0.9/0.7,(t3,t4)0.8/0.6);
```

## 17) `dnakit.microbiome.khmer_build_countgraph()`

- **Purpose:** Count the canonical k-mers of FASTA or FASTQ reads into a khmer countgraph and save it.
- **Method:** Reads the file with `khmer.ReadParser`, upper-cases each read, and feeds the reads that are at least `ksize` long and contain only A, C, G and T to `Countgraph.consume` (`SmallCountgraph` when `small_count=True`). Shorter reads and reads with other symbols are skipped and counted in the result. The table is saved with `save()` to a new path. Counts are approximate: a hash collision can only raise a count, and a k-mer and its reverse complement share one count.
- **API:** `khmer_build_countgraph(source: 'str', output_path: 'str', *, ksize: 'int' = 32, table_size: 'int' = 1000000, n_tables: 'int' = 4, small_count: 'bool' = False, max_reads: 'int' = 1000000, max_bases: 'int' = 100000000, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** needs the optional package khmer 2.1 (extra `dnakit[bio-khmer]`); no external program is needed.
- **Input:** `source` is an existing FASTA or FASTQ file; `output_path` is a new path in an existing folder; `ksize` is 1 to 32 (default 32); `table_size` is the number of buckets per table (default 1,000,000); `n_tables` is 1 to 20 (default 4); `small_count` is a boolean; `max_reads` (default 1,000,000) and `max_bases` (default 100,000,000) are explicit limits that raise `ConfigurationError` when the input exceeds them.
- **Output:** A dictionary with `status`, `source`, `output`, `kind` (`countgraph` or `smallcountgraph`), `ksize`, `table_size`, `n_tables`, `hashsizes` (the prime bucket count of each table), `occupied_bins`, `kmers_consumed`, `reads`, `used_reads`, `skipped_non_acgt_reads` and `skipped_short_reads`.
- **Example code:**

```python
import tempfile
from pathlib import Path

import khmer

from dnakit.microbiome import khmer_build_countgraph

with tempfile.TemporaryDirectory() as folder:
    output = Path(folder) / "reads.ct"
    result = khmer_build_countgraph("reads.fa", str(output), ksize=8, table_size=1000, n_tables=2, verbose=False)
    print(result["kind"], result["hashsizes"], result["occupied_bins"], result["kmers_consumed"])
    print(result["reads"], result["used_reads"], result["skipped_short_reads"], result["skipped_non_acgt_reads"])
    table = khmer.load_countgraph(str(output))
    print(table.get("ACGTACGT"), table.get("ACGTACGA"))
```

- **Example output:**

```text
countgraph [997, 991] 11 23
5 3 1 1
4 0
```

## 18) `dnakit.microbiome.khmer_build_nodegraph()`

- **Purpose:** Record which canonical k-mers occur in FASTA or FASTQ reads in a khmer nodegraph (presence table) and save it.
- **Method:** Does the same as `khmer_build_countgraph` but builds a `khmer.Nodegraph`, which stores one bit per bucket. It answers presence questions with a small false-positive rate that depends on the table size and the number of k-mers. The table is saved with `save()` to a new path and can be loaded with `khmer.load_nodegraph`.
- **API:** `khmer_build_nodegraph(source: 'str', output_path: 'str', *, ksize: 'int' = 32, table_size: 'int' = 1000000, n_tables: 'int' = 4, max_reads: 'int' = 1000000, max_bases: 'int' = 100000000, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** needs the optional package khmer 2.1 (extra `dnakit[bio-khmer]`); no external program is needed.
- **Input:** `source`, `output_path`, `ksize`, `table_size`, `n_tables`, `max_reads` and `max_bases` as in `khmer_build_countgraph` (there is no `small_count`).
- **Output:** The same dictionary as `khmer_build_countgraph`, with `kind` equal to `nodegraph`. `kmers_consumed` counts k-mers fed in, not distinct ones.
- **Example code:**

```python
import tempfile
from pathlib import Path

import khmer

from dnakit.microbiome import khmer_build_nodegraph

with tempfile.TemporaryDirectory() as folder:
    output = Path(folder) / "reads.nodegraph"
    result = khmer_build_nodegraph("reads.fa", str(output), ksize=8, table_size=1000, n_tables=2, verbose=False)
    print(result["kind"], result["hashsizes"], result["occupied_bins"], result["kmers_consumed"])
    table = khmer.load_nodegraph(str(output))
    print(table.get("ACGTACGT"), table.get("ACGTACGA"))
```

- **Example output:**

```text
nodegraph [997, 991] 11 23
1 0
```

## 19) `dnakit.microbiome.khmer_table_size_from_memory()`

- **Purpose:** Compute the buckets per table and the estimated memory of a khmer table from a memory budget or from a table size.
- **Method:** Reads khmer's `_buckets_per_byte` for the graph type (1 for countgraph, 2 for smallcountgraph, 8 for nodegraph). With `max_memory_bytes` the table size is `int(multiplier * buckets_per_byte * max_memory_bytes / n_tables)`, as in khmer's `calculate_graphsize`; with `max_tablesize` it is used as given. The memory estimate is `n_tables * table_size / buckets_per_byte`. Nothing is built or written.
- **API:** `khmer_table_size_from_memory(graph_type: "Literal['countgraph', 'smallcountgraph', 'nodegraph']", *, max_memory_bytes: 'int | None' = None, max_tablesize: 'int | None' = None, n_tables: 'int' = 4, multiplier: 'float' = 1.0, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** needs the optional package khmer 2.1 (extra `dnakit[bio-khmer]`); no external program is needed.
- **Input:** `graph_type` is `countgraph`, `smallcountgraph` or `nodegraph`; exactly one of `max_memory_bytes` (a positive number) and `max_tablesize` (an integer from 1 to 10**12) must be given; `n_tables` is 1 to 20 (default 4); `multiplier` is a positive number that scales the memory-derived size.
- **Output:** A dictionary with `status`, `operation`, `graph_type`, `buckets_per_byte`, `n_tables`, `table_size` and `estimated_memory_bytes` (a float).
- **Example code:**

```python
from dnakit.microbiome import khmer_table_size_from_memory

for graph_type in ("countgraph", "smallcountgraph", "nodegraph"):
    result = khmer_table_size_from_memory(graph_type, max_memory_bytes=1_000_000, verbose=False)
    print(result["graph_type"], result["buckets_per_byte"], result["table_size"], result["estimated_memory_bytes"])
```

- **Example output:**

```text
countgraph 1 250000 1000000.0
smallcountgraph 2 500000 1000000.0
nodegraph 8 2000000 1000000.0
```
