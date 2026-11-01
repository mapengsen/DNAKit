# gReLU: model zoo, datasets, plots and motif tools

This group of functions in `dnakit.predictions` is a set of adapters over [gReLU](https://github.com/Genentech/gReLU) (MIT licence), Genentech's library for training and interpreting sequence-to-function models. Each function calls one gReLU function with the caller's arguments and does not re-implement it. The group covers genome and interval helpers (chromosome-end checks, GC-matched background intervals, sequence extraction, variant sequences, chromosome sizes, GTF reading), the gReLU model zoo on [Hugging Face](https://huggingface.co/Genentech) and on legacy [Weights & Biases](https://wandb.ai/) projects, JASPAR motifs ([pyjaspar](https://github.com/asntech/pyjaspar)), the tomtom motif comparison ([memelite](https://github.com/jmschrei/memesuite-lite)), and plots of predictions, attributions and motifs. Genomes are resolved through [genomepy](https://github.com/vanheeringen-lab/genomepy). All of these packages are MIT or Apache-2.0 licensed.

Import the functions with `from dnakit.predictions import ...`. The dependencies are optional (extra `dnakit[bio-grelu]`: `grelu`, `huggingface-hub`, `wandb`, `memelite`) and are imported only when a function is called. If a package is missing, the function raises `BackendUnavailableError` with an install hint. Functions that reach a remote service (genome download, the Hugging Face and Weights & Biases model zoos) need internet access only when they are called. All examples on this page were run offline, with gReLU 1.1.0.

Three points apply to many functions:

- Wherever a `genome` argument is a genome name (for example `hg38`), gReLU installs it with genomepy on first use, which downloads it. If `genome` is the path of a local FASTA file instead, gReLU uses that file directly (it needs a `.fai` index, which pyfaidx builds on first use, and a `<file>.sizes` file of chromosome sizes). The examples use such a local file, so they run offline.
- gReLU's own messages (such as `Keeping 1 intervals`) are printed to the standard output. gReLU errors are re-raised as `BackendExecutionError` with the upstream message.
- Plot functions save an image (`.png`, `.pdf` or `.svg`; the file must not exist yet and its folder must exist) and return a dictionary with `status`, `backend`, `version`, `function`, `image_path` and `data` (a description of the plotted inputs). Several of gReLU's plot functions (marked in their sections) return a plotnine object that must be drawn before it is saved. DNAKit saves the current matplotlib figure without drawing it, so for those functions the saved image is blank (checked with gReLU 1.1.0 and plotnine 0.15.8).

## Example data

The examples use a two-sequence genome `tiny.fa` (40 bp and 30 bp) with its chromosome-size file `tiny.fa.sizes` (tab-separated), a small GTF file `genes.gtf` (tab-separated) and a two-motif MEME database `db.meme`. Other inputs are built in the example code.

**Example file:** `tiny.fa`

```text
>chr1
ACGTTGCAAGGCTTAACCGGTATCGATCGGATCCAGTCAA
>chr2
TTGACCAGTCAGGCATCGATTCGAGGCCTA
```

**Example file:** `tiny.fa.sizes`

```text
chr1	40
chr2	30
```

**Example file:** `genes.gtf`

```text
chr1	ex	gene	101	400	.	+	.	gene_id "g1"; gene_name "GENEA";
chr1	ex	exon	101	200	.	+	.	gene_id "g1"; gene_name "GENEA";
chr1	ex	exon	301	400	.	+	.	gene_id "g1"; gene_name "GENEA";
chr2	ex	gene	51	150	.	-	.	gene_id "g2"; gene_name "GENEB";
```

**Example file:** `db.meme`

```text
MEME version 4

ALPHABET= ACGT

strands: + -

Background letter frequencies
A 0.25 C 0.25 G 0.25 T 0.25

MOTIF GATA_like
letter-probability matrix: alength= 4 w= 4 nsites= 20 E= 0
0.05 0.05 0.85 0.05
0.85 0.05 0.05 0.05
0.05 0.05 0.05 0.85
0.85 0.05 0.05 0.05

MOTIF CACGTG_like
letter-probability matrix: alength= 4 w= 6 nsites= 20 E= 0
0.05 0.85 0.05 0.05
0.85 0.05 0.05 0.05
0.05 0.85 0.05 0.05
0.05 0.05 0.85 0.05
0.05 0.05 0.05 0.85
0.05 0.05 0.85 0.05
```

## 1) `dnakit.predictions.grelu_check_chrom_ends()`

- **Purpose:** Check that genomic intervals do not extend past the chromosome ends of a genome.
- **Method:** Calls gReLU `grelu.data.preprocess.check_chrom_ends`. It flags intervals with a negative start and, when a genome is given, intervals whose end is larger than the chromosome size from `read_sizes`. If a genome is not given, only negative starts are flagged and gReLU issues a warning. It returns nothing when all intervals are valid and raises an error listing the row labels of the failing intervals otherwise.
- **API:** `grelu_check_chrom_ends(data: 'Any', genome: 'str | None' = None, *, verbose: bool = True) -> 'None'`.
- **Requirements:** Needs the optional packages grelu and genomepy (extra `dnakit[bio-grelu]`). A genome name is downloaded on first use (needs internet access); a local FASTA path with a `.sizes` file works offline.
- **Input:** `data` is a pandas DataFrame with columns `chrom`, `start`, `end` (or an AnnData object with intervals in `.var`). `genome` is a genome name or the path of a local FASTA file.
- **Output:** `None` when every interval is inside the chromosomes. Otherwise `BackendExecutionError` whose message lists the indices of the offending intervals.
- **Example code:**

```python
import pandas as pd

from dnakit.exceptions import BackendExecutionError
from dnakit.predictions import grelu_check_chrom_ends

intervals = pd.DataFrame(
    {"chrom": ["chr1", "chr1", "chr2"], "start": [-5, 10, 20], "end": [10, 20, 35]}
)
print(grelu_check_chrom_ends(intervals.iloc[[1]], "tiny.fa", verbose=False))
try:
    grelu_check_chrom_ends(intervals, "tiny.fa", verbose=False)
except BackendExecutionError as error:
    print(type(error).__name__, str(error).split("failed: ")[1])
```

- **Example output:**

```text
None
BackendExecutionError Indices of intervals that extend beyond the chromosome ends: 0,2.
```

## 2) `dnakit.predictions.grelu_check_reference()`

- **Purpose:** Check that the reference alleles of a variant table match the reference genome.
- **Method:** Calls gReLU `grelu.variant.check_reference`. For each variant it extracts the genome base(s) at `pos` (1-based) with `intervals_to_strings` and compares them with `ref`. Rows whose `ref` equals `null_string` are skipped. A mismatch is reported as a Python `UserWarning` that lists the row labels, not as an error.
- **API:** `grelu_check_reference(variants: 'Any', *, genome: 'str' = 'hg38', null_string: 'str' = '-', verbose: bool = True) -> 'None'`.
- **Requirements:** Needs the optional packages grelu and genomepy (extra `dnakit[bio-grelu]`). A genome name is downloaded on first use; a local FASTA path works offline.
- **Input:** `variants` is a DataFrame with columns `chrom`, `pos`, `ref` and `alt`. `genome` is a genome name or a FASTA path (default `hg38`). `null_string` marks an absent base (default `-`).
- **Output:** `None`. Check the warnings (for example with `warnings.catch_warnings(record=True)`) for the mismatching rows. The input table is left unchanged.
- **Example code:**

```python
import warnings

import pandas as pd

from dnakit.predictions import grelu_check_reference

variants = pd.DataFrame(
    {"chrom": ["chr1", "chr1"], "pos": [5, 6], "ref": ["T", "A"], "alt": ["G", "C"]}
)
with warnings.catch_warnings(record=True) as caught:
    warnings.simplefilter("always")
    result = grelu_check_reference(variants, genome="tiny.fa", verbose=False)
messages = [str(w.message) for w in caught if "not centered" in str(w.message)]
print(result, len(messages))
print(messages[0])
```

- **Example output:**

```text
None 1
Sequences are not centered on reference at 1 indices [1]
```

## 3) `dnakit.predictions.grelu_download_dataset()`

- **Purpose:** Download a dataset file of the gReLU model zoo from Hugging Face and return its local path.
- **Method:** Calls gReLU `grelu.resources.download_dataset`, which calls `huggingface_hub.hf_hub_download` with `repo_type="dataset"`. The file is stored in the Hugging Face cache, and the cached path is returned.
- **API:** `grelu_download_dataset(repo_id: 'str', filename: 'str' = 'data.h5ad', *, verbose: bool = True) -> 'str'`.
- **Requirements:** Needs the optional packages grelu and huggingface_hub (extra `dnakit[bio-grelu]`) and internet access (no key for public repositories). Setting `HF_HUB_OFFLINE=1` makes it use only the local cache.
- **Input:** `repo_id` is a non-empty dataset repository ID, for example `Genentech/human-atac-catlas-data`; `filename` is the file in the repository (default `data.h5ad`).
- **Output:** The local file path as a string. Download failures are raised as `BackendExecutionError`.
- **Example code:**

```python
import os

os.environ["HF_HUB_OFFLINE"] = "1"

from dnakit.exceptions import BackendExecutionError
from dnakit.predictions import grelu_download_dataset

try:
    grelu_download_dataset("example-org/missing-data", verbose=False)
except BackendExecutionError as error:
    print(type(error).__name__, "local cache" in str(error))
```

- **Example output:**

```text
BackendExecutionError True
```

Usage with the dependency installed:

```python
from dnakit.predictions import grelu_download_dataset

path = grelu_download_dataset("Genentech/human-atac-catlas-data", "data.h5ad")
```

## 4) `dnakit.predictions.grelu_download_model()`

- **Purpose:** Download a model checkpoint of the gReLU model zoo from Hugging Face and return its local path.
- **Method:** Calls gReLU `grelu.resources.download_model`, which calls `huggingface_hub.hf_hub_download`. The file is stored in the Hugging Face cache.
- **API:** `grelu_download_model(repo_id: 'str', filename: 'str' = 'model.ckpt', *, verbose: bool = True) -> 'str'`.
- **Requirements:** Needs the optional packages grelu and huggingface_hub (extra `dnakit[bio-grelu]`) and internet access (no key for public repositories). Setting `HF_HUB_OFFLINE=1` makes it use only the local cache.
- **Input:** `repo_id` is a non-empty model repository ID, for example `Genentech/human-atac-catlas-model`; `filename` is the checkpoint file (default `model.ckpt`).
- **Output:** The local file path as a string. Failures are raised as `BackendExecutionError`.
- **Example code:**

```python
import os

os.environ["HF_HUB_OFFLINE"] = "1"

from dnakit.exceptions import BackendExecutionError
from dnakit.predictions import grelu_download_model

try:
    grelu_download_model("example-org/missing-model", verbose=False)
except BackendExecutionError as error:
    print(type(error).__name__, "local cache" in str(error))
```

- **Example output:**

```text
BackendExecutionError True
```

Usage with the dependency installed:

```python
from dnakit.predictions import grelu_download_model

path = grelu_download_model("Genentech/human-atac-catlas-model")
```

## 5) `dnakit.predictions.grelu_filter_chrom_ends()`

- **Purpose:** Drop intervals that extend beyond the chromosome ends, optionally after trimming a margin from both ends.
- **Method:** Calls gReLU `grelu.data.preprocess.filter_chrom_ends`. It keeps intervals with `start >= pad` and, when a genome is given, `end <= chromosome size - pad`. gReLU prints `Keeping N intervals`.
- **API:** `grelu_filter_chrom_ends(data: 'Any', genome: 'str | None' = None, *, pad: 'int' = 0, inplace: 'bool' = False, verbose: bool = True) -> 'Any'`.
- **Requirements:** Needs the optional packages grelu and genomepy (extra `dnakit[bio-grelu]`). A genome name is downloaded on first use; a local FASTA path with a `.sizes` file works offline.
- **Input:** `data` is a DataFrame with `chrom`, `start`, `end` (or an AnnData object with intervals in `.var`). `genome` is a genome name or FASTA path. `pad` is the number of bases to ignore at each chromosome end. With `inplace=True` the input is modified.
- **Output:** The filtered intervals in the same format (a new DataFrame unless `inplace=True`).
- **Example code:**

```python
import pandas as pd

from dnakit.predictions import grelu_filter_chrom_ends

intervals = pd.DataFrame(
    {"chrom": ["chr1", "chr1", "chr1"], "start": [3, 10, 33], "end": [10, 20, 38]}
)
print(len(grelu_filter_chrom_ends(intervals, "tiny.fa", verbose=False)))
kept = grelu_filter_chrom_ends(intervals, "tiny.fa", pad=5, verbose=False)
print(kept.to_dict("records"))
```

- **Example output:**

```text
Keeping 3 intervals
3
Keeping 1 intervals
[{'chrom': 'chr1', 'start': 10, 'end': 20}]
```

## 6) `dnakit.predictions.grelu_get_base_models()`

- **Purpose:** List the base models a model of the zoo was fine-tuned from.
- **Method:** Calls gReLU `grelu.resources.get_base_models`, which reads the `base_model` field of the model card with `HfApi().model_info`. A single string becomes a one-element list; an empty field gives an empty list.
- **API:** `grelu_get_base_models(repo_id: 'str', *, verbose: bool = True) -> 'list[str]'`.
- **Requirements:** Needs the optional packages grelu and huggingface_hub (extra `dnakit[bio-grelu]`) and internet access.
- **Input:** `repo_id` is a non-empty model repository ID.
- **Output:** A list of base-model repository IDs.
- **Example code:**

```python
import os

os.environ["HF_HUB_OFFLINE"] = "1"

from dnakit.exceptions import BackendExecutionError
from dnakit.predictions import grelu_get_base_models

try:
    grelu_get_base_models("Genentech/borzoi-model", verbose=False)
except BackendExecutionError as error:
    print(type(error).__name__, "offline mode is enabled" in str(error))
```

- **Example output:**

```text
BackendExecutionError True
```

Usage with the dependency installed:

```python
from dnakit.predictions import grelu_get_base_models

bases = grelu_get_base_models("Genentech/human-atac-catlas-model")
```

## 7) `dnakit.predictions.grelu_get_dataset_info()`

- **Purpose:** Read the dataset card metadata of a dataset repository on Hugging Face.
- **Method:** Calls gReLU `grelu.resources.get_dataset_info`, which uses `HfApi().dataset_info` and `list_repo_files(repo_type="dataset")`.
- **API:** `grelu_get_dataset_info(repo_id: 'str', *, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the optional packages grelu and huggingface_hub (extra `dnakit[bio-grelu]`) and internet access.
- **Input:** `repo_id` is a non-empty dataset repository ID.
- **Output:** A dictionary with `id`, `tags`, `card_data`, `downloads`, `last_modified` and `files`.
- **Example code:**

```python
import os

os.environ["HF_HUB_OFFLINE"] = "1"

from dnakit.exceptions import BackendExecutionError
from dnakit.predictions import grelu_get_dataset_info

try:
    grelu_get_dataset_info("Genentech/human-atac-catlas-data", verbose=False)
except BackendExecutionError as error:
    print(type(error).__name__, "offline mode is enabled" in str(error))
```

- **Example output:**

```text
BackendExecutionError True
```

Usage with the dependency installed:

```python
from dnakit.predictions import grelu_get_dataset_info

info = grelu_get_dataset_info("Genentech/human-atac-catlas-data")
print(info["files"])
```

## 8) `dnakit.predictions.grelu_get_datasets_by_model()`

- **Purpose:** List the datasets linked to a model in its model card.
- **Method:** Calls gReLU `grelu.resources.get_datasets_by_model`, which returns the `datasets` field of the model card read with `HfApi().model_info` (an empty list if the field is empty).
- **API:** `grelu_get_datasets_by_model(repo_id: 'str', *, verbose: bool = True) -> 'list[str]'`.
- **Requirements:** Needs the optional packages grelu and huggingface_hub (extra `dnakit[bio-grelu]`) and internet access.
- **Input:** `repo_id` is a non-empty model repository ID.
- **Output:** A list of dataset repository IDs.
- **Example code:**

```python
import os

os.environ["HF_HUB_OFFLINE"] = "1"

from dnakit.exceptions import BackendExecutionError
from dnakit.predictions import grelu_get_datasets_by_model

try:
    grelu_get_datasets_by_model("Genentech/human-atac-catlas-model", verbose=False)
except BackendExecutionError as error:
    print(type(error).__name__, "offline mode is enabled" in str(error))
```

- **Example output:**

```text
BackendExecutionError True
```

Usage with the dependency installed:

```python
from dnakit.predictions import grelu_get_datasets_by_model

datasets = grelu_get_datasets_by_model("Genentech/human-atac-catlas-model")
```

## 9) `dnakit.predictions.grelu_get_gc_matched_intervals()`

- **Purpose:** Sample background intervals whose GC content matches the given intervals.
- **Method:** Calls gReLU `grelu.data.preprocess.get_gc_matched_intervals`, which uses tangermeme's `extract_matching_loci` on the genome FASTA with the length of the input intervals. Afterwards gReLU removes blacklisted regions: with `blacklist=None` it looks for a bundled blacklist for the genome name and, if there is none (for example for a FASTA path), skips this step and prints a message.
- **API:** `grelu_get_gc_matched_intervals(intervals: 'Any', genome: 'str', *, binwidth: 'float' = 0.1, chroms: 'str' = 'autosomes', blacklist: 'str | None' = None, seed: 'int | None' = None, verbose: bool = True) -> 'Any'`.
- **Requirements:** Needs the optional packages grelu (with tangermeme) and genomepy (extra `dnakit[bio-grelu]`). A genome name is downloaded on first use; a local FASTA path works offline.
- **Input:** `intervals` is a DataFrame with `chrom`, `start`, `end` and one common length. `genome` is a non-empty genome name or FASTA path. `binwidth` is the GC resolution. `chroms` is `autosomes` (chr1 to chr22) or another gReLU shortcut; gReLU also accepts a list of chromosome names, which is needed for a genome that does not have chr1 to chr22. `blacklist` is a blacklist file, and `seed` fixes the sampling.
- **Output:** A DataFrame of GC-matched intervals (`chrom`, `start`, `end`), as many rows as input intervals when enough candidates exist.
- **Example code:**

```python
import contextlib
import io
import random

import pandas as pd

from dnakit.predictions import grelu_get_gc_matched_intervals

rng = random.Random(2)
with open("big.fa", "w") as handle:
    for chrom in ("chr1", "chr2"):
        sequence = "".join(rng.choice("ACGT") for _ in range(20000))
        handle.write(f">{chrom}\n")
        for i in range(0, len(sequence), 60):
            handle.write(sequence[i:i + 60] + "\n")
with open("big.fa.sizes", "w") as handle:
    handle.write("chr1\t20000\nchr2\t20000\n")

positives = pd.DataFrame(
    {"chrom": ["chr1", "chr1", "chr2"], "start": [1000, 3000, 500], "end": [1100, 3100, 600]}
)
with contextlib.redirect_stdout(io.StringIO()):
    matched = grelu_get_gc_matched_intervals(
        positives, "big.fa", chroms=["chr1", "chr2"], seed=1, verbose=False
    )
print(len(matched), list(matched.columns), (matched.end - matched.start).unique().tolist())
```

- **Example output:**

```text
3 ['chrom', 'start', 'end'] [100]
```

## 10) `dnakit.predictions.grelu_get_genome()`

- **Purpose:** Load a genome as a gReLU genome object, installing it through genomepy when a name is given.
- **Method:** Calls gReLU `grelu.io.genome.get_genome`. If the argument is an existing file path, it returns a gReLU `CustomGenome` (a pyfaidx wrapper with `get_seq`, `sizes_file` and `genome_file`). Otherwise it installs the named genome with `genomepy.install_genome` if it is not installed (this downloads it) and returns `genomepy.Genome`.
- **API:** `grelu_get_genome(genome: 'str', *, verbose: bool = True) -> 'Any'`.
- **Requirements:** Needs the optional packages grelu and genomepy (extra `dnakit[bio-grelu]`). A genome name needs internet access on first use; a local FASTA path with a `.fai` index (pyfaidx builds it) works offline.
- **Input:** `genome` is a non-empty genome name or the path of a local FASTA file.
- **Output:** The genome object. `get_seq(chrom, start, end)` takes 1-based inclusive coordinates.
- **Example code:**

```python
from dnakit.predictions import grelu_get_genome

genome = grelu_get_genome("tiny.fa", verbose=False)
print(type(genome).__name__, genome.get_seq("chr2", 1, 8))
print(genome.sizes_file.endswith("tiny.fa.sizes"))
```

- **Example output:**

```text
CustomGenome TTGACCAG
True
```

## 11) `dnakit.predictions.grelu_get_jaspar()`

- **Purpose:** Fetch motifs from the JASPAR database as position probability matrices.
- **Method:** Calls gReLU `grelu.io.motifs.get_jaspar`, which opens pyjaspar's `jaspardb(release)` and calls `fetch_motifs`. pyjaspar ships the JASPAR releases as bundled SQLite files, so no network access was needed in the test. The keys are `<matrix_id>_<name>` and the values are NumPy arrays of shape (4, length). `species` accepts a taxonomy ID, or `human` and `mouse`. Note that gReLU 1.1.0 always fetches the `CORE` collection and ignores the `collection` argument.
- **API:** `grelu_get_jaspar(release: 'str' = 'JASPAR2024', collection: 'str' = 'CORE', *, tax_group: 'str | None' = None, species: 'str | None' = None, verbose: bool = True) -> 'Any'`.
- **Requirements:** Needs the optional packages grelu and pyjaspar (a dependency of grelu; extra `dnakit[bio-grelu]`). It worked without internet access with the bundled releases.
- **Input:** `release` such as `JASPAR2024`; `tax_group` such as `plants` or `vertebrates`; `species` as a taxonomy ID string, for example `3702` for Arabidopsis thaliana.
- **Output:** A dictionary mapping motif names to arrays of shape (4, length).
- **Example code:**

```python
from dnakit.predictions import grelu_get_jaspar

motifs = grelu_get_jaspar(tax_group="plants", species="3702", verbose=False)
print(type(motifs).__name__, len(motifs))
name = sorted(motifs)[0]
print(name, motifs[name].shape)
```

- **Example output:**

```text
dict 591
MA0001.3_AGL3 (4, 10)
```

## 12) `dnakit.predictions.grelu_get_model_info()`

- **Purpose:** Read the model card metadata of a model repository on Hugging Face.
- **Method:** Calls gReLU `grelu.resources.get_model_info`, which uses `HfApi().model_info` and `list_repo_files`.
- **API:** `grelu_get_model_info(repo_id: 'str', *, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the optional packages grelu and huggingface_hub (extra `dnakit[bio-grelu]`) and internet access.
- **Input:** `repo_id` is a non-empty model repository ID.
- **Output:** A dictionary with `id`, `tags`, `card_data`, `downloads`, `last_modified` and `files`.
- **Example code:**

```python
import os

os.environ["HF_HUB_OFFLINE"] = "1"

from dnakit.exceptions import BackendExecutionError
from dnakit.predictions import grelu_get_model_info

try:
    grelu_get_model_info("Genentech/human-atac-catlas-model", verbose=False)
except BackendExecutionError as error:
    print(type(error).__name__, "offline mode is enabled" in str(error))
```

- **Example output:**

```text
BackendExecutionError True
```

Usage with the dependency installed:

```python
from dnakit.predictions import grelu_get_model_info

info = grelu_get_model_info("Genentech/human-atac-catlas-model")
print(info["files"])
```

## 13) `dnakit.predictions.grelu_get_models_by_dataset()`

- **Purpose:** List the models of the gReLU collection that were trained on a dataset.
- **Method:** Calls gReLU `grelu.resources.get_models_by_dataset`. It lists the models of the model-zoo collection and, for each one, reads the `datasets` field of its model card, so it sends one request per model.
- **API:** `grelu_get_models_by_dataset(repo_id: 'str', *, verbose: bool = True) -> 'list[str]'`.
- **Requirements:** Needs the optional packages grelu and huggingface_hub (extra `dnakit[bio-grelu]`) and internet access.
- **Input:** `repo_id` is a non-empty dataset repository ID.
- **Output:** A list of model repository IDs.
- **Example code:**

```python
import os

os.environ["HF_HUB_OFFLINE"] = "1"

from dnakit.exceptions import BackendExecutionError
from dnakit.predictions import grelu_get_models_by_dataset

try:
    grelu_get_models_by_dataset("Genentech/human-atac-catlas-data", verbose=False)
except BackendExecutionError as error:
    print(type(error).__name__, "offline mode is enabled" in str(error))
```

- **Example output:**

```text
BackendExecutionError True
```

Usage with the dependency installed:

```python
from dnakit.predictions import grelu_get_models_by_dataset

models = grelu_get_models_by_dataset("Genentech/human-atac-catlas-data")
```

## 14) `dnakit.predictions.grelu_intervals_to_strings()`

- **Purpose:** Extract the DNA sequences of genomic intervals from a genome.
- **Method:** Calls gReLU `grelu.sequence.format.intervals_to_strings`. Intervals are 0-based half-open (`start`, `end`); sequences are upper-cased and reverse-complemented for rows whose `strand` is `-` when a `strand` column exists.
- **API:** `grelu_intervals_to_strings(intervals: 'Any', genome: 'str', *, verbose: bool = True) -> 'Any'`.
- **Requirements:** Needs the optional packages grelu and genomepy (extra `dnakit[bio-grelu]`). A genome name is downloaded on first use; a local FASTA path works offline.
- **Input:** `intervals` is a DataFrame with `chrom`, `start`, `end` (and optionally `strand`), a Series or a dictionary for one interval. `genome` is a non-empty genome name or FASTA path.
- **Output:** A list of sequence strings for a DataFrame, or one string for a Series or dictionary.
- **Example code:**

```python
import pandas as pd

from dnakit.predictions import grelu_intervals_to_strings

intervals = pd.DataFrame(
    {"chrom": ["chr1", "chr2"], "start": [10, 20], "end": [20, 30], "strand": ["+", "-"]}
)
print(grelu_intervals_to_strings(intervals, "tiny.fa", verbose=False))
print(grelu_intervals_to_strings(intervals.iloc[0], "tiny.fa", verbose=False))
```

- **Example output:**

```text
['GCTTAACCGG', 'TAGGCCTCGA']
GCTTAACCGG
```

## 15) `dnakit.predictions.grelu_list_datasets()`

- **Purpose:** List the dataset repository IDs of the gReLU model zoo on Hugging Face.
- **Method:** Calls gReLU `grelu.resources.list_datasets`, which reads the model-zoo collection with `HfApi().get_collection` and keeps the items whose ID ends with `-data`.
- **API:** `grelu_list_datasets(*, verbose: bool = True) -> 'list[str]'`.
- **Requirements:** Needs the optional packages grelu and huggingface_hub (extra `dnakit[bio-grelu]`) and internet access.
- **Input:** None.
- **Output:** A list of dataset repository IDs.
- **Example code:**

```python
import os

os.environ["HF_HUB_OFFLINE"] = "1"

from dnakit.exceptions import BackendExecutionError
from dnakit.predictions import grelu_list_datasets

try:
    grelu_list_datasets(verbose=False)
except BackendExecutionError as error:
    print(type(error).__name__, "offline mode is enabled" in str(error))
```

- **Example output:**

```text
BackendExecutionError True
```

Usage with the dependency installed:

```python
from dnakit.predictions import grelu_list_datasets

datasets = grelu_list_datasets()
```

## 16) `dnakit.predictions.grelu_list_models()`

- **Purpose:** List the model repository IDs of the gReLU model zoo on Hugging Face.
- **Method:** Calls gReLU `grelu.resources.list_models`, which reads the model-zoo collection with `HfApi().get_collection` and keeps the items whose ID ends with `-model`.
- **API:** `grelu_list_models(*, verbose: bool = True) -> 'list[str]'`.
- **Requirements:** Needs the optional packages grelu and huggingface_hub (extra `dnakit[bio-grelu]`) and internet access.
- **Input:** None.
- **Output:** A list of model repository IDs.
- **Example code:**

```python
import os

os.environ["HF_HUB_OFFLINE"] = "1"

from dnakit.exceptions import BackendExecutionError
from dnakit.predictions import grelu_list_models

try:
    grelu_list_models(verbose=False)
except BackendExecutionError as error:
    print(type(error).__name__, "offline mode is enabled" in str(error))
```

- **Example output:**

```text
BackendExecutionError True
```

Usage with the dependency installed:

```python
from dnakit.predictions import grelu_list_models

models = grelu_list_models()
```

## 17) `dnakit.predictions.grelu_load_model()`

- **Purpose:** Load a gReLU model from the Hugging Face model zoo.
- **Method:** Calls gReLU `grelu.resources.load_model`, which downloads the checkpoint with `download_model` and loads it with `LightningModel.load_from_checkpoint(path, map_location=device)`. The `project` and `model_name` arguments are only kept to give a clear error: in gReLU 1.1.0, passing either one raises a deprecation error that points to `grelu_wandb_load_model`. The DNAKit docstring says that this function can load from Weights & Biases, but it cannot.
- **API:** `grelu_load_model(repo_id: 'str | None' = None, filename: 'str' = 'model.ckpt', *, device: 'str | int' = 'cpu', project: 'str | None' = None, model_name: 'str | None' = None, verbose: bool = True) -> 'Any'`.
- **Requirements:** Needs the optional packages grelu and huggingface_hub (extra `dnakit[bio-grelu]`), internet access for the first download and a model checkpoint (downloaded from the Hub).
- **Input:** `repo_id` is a model repository ID (required in practice); `filename` is the checkpoint (default `model.ckpt`); `device` is `cpu` or a GPU index.
- **Output:** A gReLU `LightningModel`. A missing `repo_id` or the legacy arguments raise `BackendExecutionError`.
- **Example code:**

```python
from dnakit.exceptions import BackendExecutionError
from dnakit.predictions import grelu_load_model

try:
    grelu_load_model(project="zoo", model_name="model", verbose=False)
except BackendExecutionError as error:
    print(type(error).__name__, "API has changed" in str(error))
try:
    grelu_load_model(verbose=False)
except BackendExecutionError as error:
    print(type(error).__name__, "repo_id is required" in str(error))
```

- **Example output:**

```text
BackendExecutionError True
BackendExecutionError True
```

Usage with the dependency installed:

```python
from dnakit.predictions import grelu_load_model

model = grelu_load_model("Genentech/human-atac-catlas-model", device="cpu")
```

## 18) `dnakit.predictions.grelu_make_insertion_bigwig()`

- **Purpose:** Convert a fragment file of Tn5 insertions into a strand-specific insertion bigWig.
- **Method:** Calls gReLU `grelu.data.preprocess.make_insertion_bigwig`. It runs a shell pipeline (`cat` or `zcat`, `awk`, `sort`, optional `grep` for `chroms`, `bedtools genomecov -bg -5`, `bedtools sort`) to make a bedGraph in `tmp_dir`, then `bedGraphToBigWig` with the genome's chromosome-size file to make the bigWig in `out_dir`, and deletes the bedGraph. DNAKit first checks that `frag_file` exists.
- **API:** `grelu_make_insertion_bigwig(frag_file: 'str', genome: 'str', *, out_prefix: 'str | None' = None, plus_shift: 'int' = 0, minus_shift: 'int' = 0, chroms: 'Sequence[str] | str | None' = None, tmp_dir: 'str' = './', out_dir: 'str' = './', verbose: bool = True) -> 'str'`.
- **Requirements:** Needs the optional packages grelu and genomepy (extra `dnakit[bio-grelu]`) and the external programs `bedtools`, `awk`, `sort`, `cat`/`zcat` and the UCSC tool `bedGraphToBigWig` on PATH. `bedGraphToBigWig` was not installed in the test environment, so the example shows the input check; the genome name is downloaded on first use.
- **Input:** `frag_file` is an existing fragment file (BED-like, optionally `.gz`); `genome` is a genome name or FASTA path; `out_prefix` is the output file name (the default is the fragment file name plus `.bw`; despite its name it is used as the full file name); `plus_shift` and `minus_shift` shift the insertion sites; `chroms` limits the chromosomes; `tmp_dir` and `out_dir` are existing folders.
- **Output:** The path of the bigWig file as a string.
- **Example code:**

```python
from dnakit.exceptions import ConfigurationError
from dnakit.predictions import grelu_make_insertion_bigwig

try:
    grelu_make_insertion_bigwig("missing_fragments.tsv", "tiny.fa", verbose=False)
except ConfigurationError as error:
    print(type(error).__name__, str(error).split(":")[0])
```

- **Example output:**

```text
ConfigurationError [CONFIGURATION_ERROR] frag_file must be an existing file
```

Usage with the dependency installed:

```python
from dnakit.predictions import grelu_make_insertion_bigwig

bigwig = grelu_make_insertion_bigwig(
    "fragments.tsv.gz", "hg38", out_prefix="insertions.bw", plus_shift=4, minus_shift=-5
)
```

## 19) `dnakit.predictions.grelu_plot_attention_matrix()`

- **Purpose:** Draw a heatmap of a transformer attention matrix over a sequence window.
- **Method:** Calls gReLU `grelu.visualize.plot_attention_matrix`, which builds a seaborn heatmap with the bin coordinates between `start_pos` and `end_pos` and optionally draws a square for each highlight interval. The figure is saved to `output_path`.
- **API:** `grelu_plot_attention_matrix(attn: 'Any', output_path: 'str', *, start_pos: 'int' = 0, end_pos: 'int | None' = None, highlight_intervals: 'Any' = None, figsize: 'tuple[float, float]' = (5, 4), verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the optional package grelu (extra `dnakit[bio-grelu]`) and matplotlib.
- **Input:** `attn` is a finite square 2-D array; `output_path` is a new `.png`, `.pdf` or `.svg` file in an existing folder; `start_pos` and `end_pos` give the genomic window (bin size is `(end_pos - start_pos) / n`); `highlight_intervals` is a DataFrame with `start` and `end`.
- **Output:** A dictionary with `status`, `backend`, `version`, `function`, `image_path` and `data` (`{"attn": {"shape": ..., "values": ...}}`).
- **Example code:**

```python
import tempfile
from pathlib import Path

import numpy as np

from dnakit.predictions import grelu_plot_attention_matrix

attn = np.eye(6) * 0.5 + 0.1
with tempfile.TemporaryDirectory() as folder:
    result = grelu_plot_attention_matrix(attn, str(Path(folder) / "attention.png"), verbose=False)
    print(result["status"], result["function"], Path(result["image_path"]).is_file())
print(result["data"]["attn"]["shape"])
```

- **Example output:**

```text
completed grelu.visualize.plot_attention_matrix True
[6, 6]
```

## 20) `dnakit.predictions.grelu_plot_attributions()`

- **Purpose:** Plot per-base attribution scores of one sequence as a sequence logo with optional highlights.
- **Method:** Calls gReLU `grelu.visualize.plot_attributions`, which slices the (4, L) array to `start_pos:end_pos`, and draws it with logomaker as a logo in which the height of each letter is its score. Highlight boxes are added with `highlight_centers` and `highlight_width`, or `highlight_positions`. The default `end_pos=-1` drops the last position.
- **API:** `grelu_plot_attributions(attrs: 'Any', output_path: 'str', *, start_pos: 'int' = 0, end_pos: 'int' = -1, figsize: 'tuple[float, float]' = (20, 2), ticks: 'int' = 10, highlight_centers: 'Sequence[int] | None' = None, highlight_width: 'int | Sequence[int]' = 5, highlight_positions: 'Sequence[int] | None' = None, ylim: 'tuple[float, float] | None' = None, facecolor: 'str | None' = 'yellow', edgecolor: 'str | None' = None, alpha: 'float | None' = 0.15, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the optional package grelu (extra `dnakit[bio-grelu]`), with logomaker and matplotlib.
- **Input:** `attrs` is a finite 2-D array of shape (4, L), rows A, C, G, T. `output_path` is a new image file. The other arguments set the window, ticks, highlight positions, colours and axis limits.
- **Output:** The plot dictionary (see above); `data["attrs"]` holds the shape and values.
- **Example code:**

```python
import tempfile
from pathlib import Path

import numpy as np

from dnakit.predictions import grelu_plot_attributions

attrs = np.zeros((4, 12))
attrs[2, 3:7] = [0.2, 0.9, 0.5, 0.7]
attrs[0, 7:10] = [0.4, 0.6, 0.3]
with tempfile.TemporaryDirectory() as folder:
    result = grelu_plot_attributions(
        attrs, str(Path(folder) / "attributions.png"), end_pos=12,
        highlight_centers=[5], figsize=(6, 1.5), verbose=False,
    )
    print(result["status"], Path(result["image_path"]).is_file())
print(result["data"]["attrs"]["shape"])
```

- **Example output:**

```text
completed True
[4, 12]
```

## 21) `dnakit.predictions.grelu_plot_binary_preds()`

- **Purpose:** Plot predicted probabilities of binary labels as box plots per class and task.
- **Method:** Calls gReLU `grelu.visualize.plot_binary_preds`, which flattens predictions and labels to one row per observation and position, and builds a plotnine box plot of the prediction by label for each task. Known issue: gReLU returns the plotnine object without drawing it, so DNAKit saves a blank image (checked with gReLU 1.1.0).
- **API:** `grelu_plot_binary_preds(preds: 'Any', labels: 'Any', output_path: 'str', *, tasks: 'Sequence[int] | Sequence[str] | None' = None, bins: 'Sequence[int] | None' = None, figsize: 'tuple[float, float]' = (4, 3), verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the optional package grelu (extra `dnakit[bio-grelu]`), with plotnine and matplotlib.
- **Input:** `preds` and `labels` are finite arrays of the same shape (observations, tasks, length). `tasks` selects task indices or names; `bins` selects length positions; `output_path` is a new image file.
- **Output:** The plot dictionary; `data` holds `preds` and `labels` descriptions.
- **Example code:**

```python
import tempfile
from pathlib import Path

import numpy as np

from dnakit.predictions import grelu_plot_binary_preds

rng = np.random.default_rng(0)
preds = rng.random((8, 2, 3))
labels = (rng.random((8, 2, 3)) > 0.5).astype(float)
with tempfile.TemporaryDirectory() as folder:
    result = grelu_plot_binary_preds(preds, labels, str(Path(folder) / "binary.png"), verbose=False)
    print(result["status"], Path(result["image_path"]).is_file())
print(result["data"]["preds"]["shape"], result["data"]["labels"]["shape"])
```

- **Example output:**

```text
completed True
[8, 2, 3] [8, 2, 3]
```

## 22) `dnakit.predictions.grelu_plot_calibration_curve()`

- **Purpose:** Plot the calibration curve of predicted probabilities against observed binary labels.
- **Method:** Calls gReLU `grelu.visualize.plot_calibration_curve`, which bins the probabilities and plots the observed fraction of positives per bin (aggregated over tasks when `aggregate=True`). Known issue: the plotnine object is not drawn before saving, so the saved image is blank (checked with gReLU 1.1.0).
- **API:** `grelu_plot_calibration_curve(probs: 'Any', labels: 'Any', output_path: 'str', *, tasks: 'Sequence[int] | Sequence[str] | None' = None, bins: 'Sequence[int] | None' = None, aggregate: 'bool' = True, figsize: 'tuple[float, float]' = (4, 3), show_legend: 'bool' = True, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the optional package grelu (extra `dnakit[bio-grelu]`), with plotnine and matplotlib.
- **Input:** `probs` and `labels` are finite 3-D arrays of the same shape (observations, tasks, length). DNAKit also accepts a 2-D `probs`, but gReLU 1.1.0 then fails with `BackendExecutionError`. `output_path` is a new image file.
- **Output:** The plot dictionary; `data` holds `probs` and `labels` descriptions.
- **Example code:**

```python
import tempfile
from pathlib import Path

import numpy as np

from dnakit.predictions import grelu_plot_calibration_curve

rng = np.random.default_rng(1)
probs = rng.random((20, 2, 3))
labels = (rng.random((20, 2, 3)) < probs).astype(float)
with tempfile.TemporaryDirectory() as folder:
    result = grelu_plot_calibration_curve(
        probs, labels, str(Path(folder) / "calibration.png"), verbose=False
    )
    print(result["status"], Path(result["image_path"]).is_file())
print(result["data"]["probs"]["shape"])
```

- **Example output:**

```text
completed True
[20, 2, 3]
```

## 23) `dnakit.predictions.grelu_plot_distribution()`

- **Purpose:** Plot the distribution of a list of values as a histogram or density.
- **Method:** Calls gReLU `grelu.visualize.plot_distribution`, which builds a plotnine histogram or density plot. The docstring of the DNAKit function also names a box plot, but gReLU supports only `histogram` and `density`; any other `method` produces an empty plot. Known issue: the plotnine object is not drawn before saving, so the saved image is blank (checked with gReLU 1.1.0).
- **API:** `grelu_plot_distribution(values: 'Sequence[float]', output_path: 'str', *, title: 'str' = 'metric', method: 'str' = 'histogram', figsize: 'tuple[float, float]' = (4, 3), verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the optional package grelu (extra `dnakit[bio-grelu]`), with plotnine and matplotlib.
- **Input:** `values` is a 1-D sequence of numbers; `title` is the axis label; `method` is `histogram` or `density`; `output_path` is a new image file.
- **Output:** The plot dictionary; `data["values"]` is the list of values.
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.predictions import grelu_plot_distribution

values = [0.1, 0.4, 0.35, 0.8, 0.6, 0.5, 0.45, 0.9]
with tempfile.TemporaryDirectory() as folder:
    result = grelu_plot_distribution(
        values, str(Path(folder) / "distribution.png"), title="auroc", method="density",
        verbose=False,
    )
    print(result["status"], Path(result["image_path"]).is_file())
print(len(result["data"]["values"]))
```

- **Example output:**

```text
completed True
8
```

## 24) `dnakit.predictions.grelu_plot_evolution()`

- **Purpose:** Plot the change in scores over the rounds of a directed-evolution run.
- **Method:** Calls gReLU `grelu.visualize.plot_evolution`, which takes the columns `iter`, `prediction_score`, `seq_score`, `total_score` and every column after the ninth, and draws one box plot per iteration for each score (the prediction and sequence scores are dropped when `seq_score` is all zero). Known issue: the plotnine object is not drawn before saving, so the saved image is blank (checked with gReLU 1.1.0).
- **API:** `grelu_plot_evolution(df: 'Any', output_path: 'str', *, figsize: 'tuple[float, float]' = (4, 3), verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the optional package grelu (extra `dnakit[bio-grelu]`), with plotnine and matplotlib.
- **Input:** `df` is a DataFrame in the format produced by `grelu.design.evolve`, with the four score columns above and at least nine columns in total. `output_path` is a new image file.
- **Output:** The plot dictionary; `data["df"]` has `columns`, `rows` and `values`.
- **Example code:**

```python
import tempfile
from pathlib import Path

import pandas as pd

from dnakit.predictions import grelu_plot_evolution

df = pd.DataFrame(
    {
        "iter": [0, 0, 1, 1],
        "prediction_score": [0.1, 0.2, 0.4, 0.5],
        "seq_score": [0.0, 0.0, 0.0, 0.0],
        "total_score": [0.1, 0.2, 0.4, 0.5],
        "seq": ["ACGT"] * 4,
        "chrom": ["x"] * 4,
        "start": [0] * 4,
        "end": [4] * 4,
        "strand": ["+"] * 4,
    }
)
with tempfile.TemporaryDirectory() as folder:
    result = grelu_plot_evolution(df, str(Path(folder) / "evolution.png"), verbose=False)
    print(result["status"], Path(result["image_path"]).is_file())
print(result["data"]["df"]["rows"], result["data"]["df"]["columns"][:2])
```

- **Example output:**

```text
completed True
4 ['iter', 'prediction_score']
```

## 25) `dnakit.predictions.grelu_plot_gc_match()`

- **Purpose:** Plot histograms of the GC content of positive and GC-matched negative intervals.
- **Method:** Calls gReLU `grelu.visualize.plot_gc_match`, which computes the GC distribution of each interval set with `gc_distribution` on the genome and draws a dodged bar plot with plotnine. Known issue: the plotnine object is not drawn before saving, so the saved image is blank (checked with gReLU 1.1.0).
- **API:** `grelu_plot_gc_match(positives: 'Any', negatives: 'Any', output_path: 'str', *, binwidth: 'float' = 0.1, genome: 'str' = 'hg38', figsize: 'tuple[float, float]' = (4, 3), verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the optional packages grelu, plotnine and genomepy (extra `dnakit[bio-grelu]`). A genome name is downloaded on first use; a local FASTA path works offline.
- **Input:** `positives` and `negatives` are DataFrames with `chrom`, `start`, `end`; `binwidth` is the GC resolution; `genome` is a genome name or FASTA path (default `hg38`); `output_path` is a new image file.
- **Output:** The plot dictionary; `data` describes the two interval tables (`columns`, `rows`, `values`).
- **Example code:**

```python
import tempfile
from pathlib import Path

import pandas as pd

from dnakit.predictions import grelu_plot_gc_match

positives = pd.DataFrame({"chrom": ["chr1", "chr1"], "start": [0, 10], "end": [10, 20]})
negatives = pd.DataFrame({"chrom": ["chr1", "chr2"], "start": [20, 0], "end": [30, 10]})
with tempfile.TemporaryDirectory() as folder:
    result = grelu_plot_gc_match(
        positives, negatives, str(Path(folder) / "gc.png"), binwidth=0.25,
        genome="tiny.fa", verbose=False,
    )
    print(result["status"], Path(result["image_path"]).is_file())
print(result["data"]["positives"]["rows"], result["data"]["negatives"]["rows"])
```

- **Example output:**

```text
completed True
2 2
```

## 26) `dnakit.predictions.grelu_plot_ism()`

- **Purpose:** Plot an in-silico mutagenesis table as a heatmap or a sequence logo.
- **Method:** Calls gReLU `grelu.visualize.plot_ISM`. With `method="heatmap"` it draws a seaborn heatmap (rows A, C, G, T; one column per position) of the table region `start_pos:end_pos`. With `method="logo"` it draws a logo of the negative mean mutation effect at each position, placed on the reference base given by the column name. The DNAKit docstring says "line plot", but gReLU supports only `heatmap` and `logo`; any other value fails with `BackendExecutionError`.
- **API:** `grelu_plot_ism(ism_preds: 'Any', output_path: 'str', *, start_pos: 'int | None' = None, end_pos: 'int | None' = None, figsize: 'tuple[float, float]' = (8, 1.5), method: 'str' = 'heatmap', verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the optional package grelu (extra `dnakit[bio-grelu]`), with seaborn, logomaker and matplotlib.
- **Input:** `ism_preds` is a DataFrame with the bases A, C, G, T as the row index and one column per sequence position, named with the reference base (as produced by `grelu.model.interpret.ISM_predict`). `start_pos` and `end_pos` select columns; `method` is `heatmap` or `logo`; `output_path` is a new image file.
- **Output:** The plot dictionary; `data["ism_preds"]` has `columns`, `rows` and `values`.
- **Example code:**

```python
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd

from dnakit.predictions import grelu_plot_ism

ism = pd.DataFrame(
    np.arange(24).reshape(4, 6) / 10, index=list("ACGT"), columns=list("ACGTAC")
)
with tempfile.TemporaryDirectory() as folder:
    heat = grelu_plot_ism(ism, str(Path(folder) / "ism.png"), verbose=False)
    logo = grelu_plot_ism(ism, str(Path(folder) / "ism_logo.svg"), method="logo", verbose=False)
    print(heat["status"], logo["status"], Path(logo["image_path"]).is_file())
print(heat["data"]["ism_preds"]["rows"], heat["data"]["ism_preds"]["columns"])
```

- **Example output:**

```text
completed completed True
4 ['A', 'C', 'G', 'T', 'A', 'C']
```

## 27) `dnakit.predictions.grelu_plot_motif()`

- **Purpose:** Draw a sequence logo of a motif given as a position matrix with four columns.
- **Method:** Calls gReLU `grelu.visualize.plot_motif`, which passes the transposed matrix to tangermeme's `plot_pwm`. The figure is saved to `output_path`.
- **API:** `grelu_plot_motif(motif: 'Any', output_path: 'str', *, name: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the optional package grelu (extra `dnakit[bio-grelu]`), with tangermeme and matplotlib.
- **Input:** `motif` is a finite 2-D array of shape (length, 4), columns A, C, G, T (probabilities or scores). `name` is an optional title. `output_path` is a new image file.
- **Output:** The plot dictionary; `data["motif"]` holds the shape and values.
- **Example code:**

```python
import tempfile
from pathlib import Path

import numpy as np

from dnakit.predictions import grelu_plot_motif

motif = np.array([[0.7, 0.1, 0.1, 0.1], [0.1, 0.7, 0.1, 0.1], [0.1, 0.1, 0.7, 0.1]])
with tempfile.TemporaryDirectory() as folder:
    result = grelu_plot_motif(motif, str(Path(folder) / "motif.svg"), name="demo", verbose=False)
    print(result["status"], Path(result["image_path"]).is_file())
print(result["data"]["motif"]["shape"])
```

- **Example output:**

```text
completed True
[3, 4]
```

## 28) `dnakit.predictions.grelu_plot_position_effect()`

- **Purpose:** Plot how the model output changes with the position of an inserted motif.
- **Method:** Calls gReLU `grelu.visualize.plot_position_effect`, which computes the mean, the standard error and a 95 % confidence interval per position, and draws them with plotnine. Known issue: the plotnine object is not drawn before saving, so the saved image is blank (checked with gReLU 1.1.0).
- **API:** `grelu_plot_position_effect(preds: 'Any', positions: 'Sequence[int]', output_path: 'str', *, title: 'str | None' = None, xlab: 'str | None' = None, figsize: 'tuple[float, float]' = (6, 3), verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the optional package grelu (extra `dnakit[bio-grelu]`), with plotnine, scipy and matplotlib.
- **Input:** `preds` is a finite 2-D array (sequences, positions), as produced by `marginalize_pattern_spacing` or `shuffle_tiles`; `positions` is a list with one value per column; `title` and `xlab` are optional labels; `output_path` is a new image file.
- **Output:** The plot dictionary; `data` holds `preds` and `positions`.
- **Example code:**

```python
import tempfile
from pathlib import Path

import numpy as np

from dnakit.predictions import grelu_plot_position_effect

preds = np.random.default_rng(0).random((5, 6))
with tempfile.TemporaryDirectory() as folder:
    result = grelu_plot_position_effect(
        preds, [0, 10, 20, 30, 40, 50], str(Path(folder) / "position.png"),
        title="motif effect", xlab="position", verbose=False,
    )
    print(result["status"], Path(result["image_path"]).is_file())
print(result["data"]["preds"]["shape"], result["data"]["positions"])
```

- **Example output:**

```text
completed True
[5, 6] [0, 10, 20, 30, 40, 50]
```

## 29) `dnakit.predictions.grelu_plot_pred_distribution()`

- **Purpose:** Plot the distribution of model predictions for the positive and negative labels of each task.
- **Method:** Calls gReLU `grelu.visualize.plot_pred_distribution`, which collects predictions and labels per task and draws density plots with plotnine. Known issue: the plotnine object is not drawn before saving, so the saved image is blank (checked with gReLU 1.1.0).
- **API:** `grelu_plot_pred_distribution(preds: 'Any', labels: 'Any', output_path: 'str', *, tasks: 'Sequence[int] | Sequence[str] | None' = None, bins: 'Sequence[int] | None' = None, figsize: 'tuple[float, float]' = (4, 3), verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the optional package grelu (extra `dnakit[bio-grelu]`), with plotnine and matplotlib.
- **Input:** `preds` and `labels` are finite 3-D arrays of the same shape (observations, tasks, length); `tasks` and `bins` select tasks and positions; `output_path` is a new image file.
- **Output:** The plot dictionary; `data` holds `preds` and `labels` descriptions.
- **Example code:**

```python
import tempfile
from pathlib import Path

import numpy as np

from dnakit.predictions import grelu_plot_pred_distribution

rng = np.random.default_rng(2)
preds = rng.random((10, 2, 3))
labels = (rng.random((10, 2, 3)) > 0.5).astype(float)
with tempfile.TemporaryDirectory() as folder:
    result = grelu_plot_pred_distribution(
        preds, labels, str(Path(folder) / "pred_distribution.png"), tasks=[0], verbose=False
    )
    print(result["status"], Path(result["image_path"]).is_file())
print(result["data"]["labels"]["shape"])
```

- **Example output:**

```text
completed True
[10, 2, 3]
```

## 30) `dnakit.predictions.grelu_plot_pred_scatter()`

- **Purpose:** Plot predictions against measured labels for each task as a scatter plot.
- **Method:** Calls gReLU `grelu.visualize.plot_pred_scatter`, which collects predictions and labels per task and draws a scatter plot with plotnine (density-coloured when `density=True`). Known issue: the plotnine object is not drawn before saving, so the saved image is blank (checked with gReLU 1.1.0).
- **API:** `grelu_plot_pred_scatter(preds: 'Any', labels: 'Any', output_path: 'str', *, tasks: 'Sequence[int] | Sequence[str] | None' = None, bins: 'Sequence[int] | None' = None, density: 'bool' = False, figsize: 'tuple[float, float]' = (4, 3), verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the optional package grelu (extra `dnakit[bio-grelu]`), with plotnine and matplotlib.
- **Input:** `preds` and `labels` are finite 3-D arrays of the same shape (observations, tasks, length); `tasks`, `bins` and `density` control what is drawn; `output_path` is a new image file.
- **Output:** The plot dictionary; `data` holds `preds` and `labels` descriptions.
- **Example code:**

```python
import tempfile
from pathlib import Path

import numpy as np

from dnakit.predictions import grelu_plot_pred_scatter

rng = np.random.default_rng(3)
labels = rng.random((12, 2, 1))
preds = labels + rng.normal(scale=0.1, size=labels.shape)
with tempfile.TemporaryDirectory() as folder:
    result = grelu_plot_pred_scatter(preds, labels, str(Path(folder) / "scatter.png"), verbose=False)
    print(result["status"], Path(result["image_path"]).is_file())
print(result["data"]["preds"]["shape"])
```

- **Example output:**

```text
completed True
[12, 2, 1]
```

## 31) `dnakit.predictions.grelu_plot_tracks()`

- **Purpose:** Plot one or more prediction or signal tracks over a sequence window.
- **Method:** Calls gReLU `grelu.visualize.plot_tracks`, which draws one stacked panel per track with matplotlib between `start_pos` and `end_pos`, optional panel titles and optional highlight intervals.
- **API:** `grelu_plot_tracks(tracks: 'Any', output_path: 'str', *, start_pos: 'int' = 0, end_pos: 'int | None' = None, titles: 'Sequence[str] | None' = None, figsize: 'tuple[float, float]' = (20, 1.5), highlight_intervals: 'Any' = None, facecolor: 'str | None' = 'yellow', edgecolor: 'str | None' = None, alpha: 'float | None' = 0.15, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the optional package grelu (extra `dnakit[bio-grelu]`), with matplotlib.
- **Input:** `tracks` is a finite 2-D array (tracks, positions); `titles` has one label per track; `highlight_intervals` is a DataFrame with `start` and `end`; `output_path` is a new image file.
- **Output:** The plot dictionary; `data["tracks"]` holds the shape and values.
- **Example code:**

```python
import tempfile
from pathlib import Path

import numpy as np

from dnakit.predictions import grelu_plot_tracks

x = np.linspace(0, 6, 40)
tracks = np.vstack([np.sin(x) + 1, np.cos(x) + 1])
with tempfile.TemporaryDirectory() as folder:
    result = grelu_plot_tracks(
        tracks, str(Path(folder) / "tracks.png"), titles=["sin", "cos"], figsize=(6, 3),
        verbose=False,
    )
    print(result["status"], Path(result["image_path"]).is_file())
print(result["data"]["tracks"]["shape"])
```

- **Example output:**

```text
completed True
[2, 40]
```

## 32) `dnakit.predictions.grelu_read_gtf()`

- **Purpose:** Read the gene features of a genome annotation as a table.
- **Method:** Calls gReLU `grelu.io.genome.read_gtf`, which opens `genomepy.Annotation(genome)` and reads its `named_gtf` table. If the annotation of a genome name is not installed, gReLU installs it with genomepy (a download, which may need UCSC tools). genomepy also accepts the path of a GTF file, which is read directly; this was used in the example. The table is reset to a default index and `seqname` is renamed to `chrom`, then `chrom`, `start`, `end` are moved first. `features` filters the `feature` column.
- **API:** `grelu_read_gtf(genome: 'str', features: 'Sequence[str] | str | None' = None, *, verbose: bool = True) -> 'Any'`.
- **Requirements:** Needs the optional packages grelu and genomepy (extra `dnakit[bio-grelu]`). A genome name needs internet access on first use; a local GTF path works offline.
- **Input:** `genome` is a non-empty genome name or the path of a GTF file; `features` is a feature name or list, such as `exon`, `CDS` or `transcript`.
- **Output:** A pandas DataFrame with the columns `chrom`, `start`, `end`, the name column, `source`, `feature`, `score`, `strand`, `frame` and `attribute`. The genomepy log messages about missing sibling files go to the log, not to the standard output.
- **Example code:**

```python
from dnakit.predictions import grelu_read_gtf

table = grelu_read_gtf("genes.gtf", verbose=False)
print(table.shape, table.columns.tolist()[:3])
exons = grelu_read_gtf("genes.gtf", "exon", verbose=False)
print(exons[["chrom", "start", "end"]].to_dict("records"))
```

- **Example output:**

```text
(4, 10) ['chrom', 'start', 'end']
[{'chrom': 'chr1', 'start': 101, 'end': 200}, {'chrom': 'chr1', 'start': 301, 'end': 400}]
```

## 33) `dnakit.predictions.grelu_read_sizes()`

- **Purpose:** Read the chromosome sizes of a genome as a table.
- **Method:** Calls gReLU `grelu.io.genome.read_sizes`, which loads the genome with `get_genome` (installing it with genomepy if a name is given) and reads its `sizes_file`.
- **API:** `grelu_read_sizes(genome: 'str' = 'hg38', *, verbose: bool = True) -> 'Any'`.
- **Requirements:** Needs the optional packages grelu and genomepy (extra `dnakit[bio-grelu]`). A genome name needs internet access on first use; a local FASTA path with a `.sizes` file works offline.
- **Input:** `genome` is a non-empty genome name (default `hg38`) or a FASTA path.
- **Output:** A DataFrame with the columns `chrom` and `size`.
- **Example code:**

```python
from dnakit.predictions import grelu_read_sizes

sizes = grelu_read_sizes("tiny.fa", verbose=False)
print(sizes.to_dict("records"))
```

- **Example output:**

```text
[{'chrom': 'chr1', 'size': 40}, {'chrom': 'chr2', 'size': 30}]
```

## 34) `dnakit.predictions.grelu_run_tomtom()`

- **Purpose:** Compare query motifs with a reference MEME database using the tomtom algorithm.
- **Method:** Calls gReLU `grelu.interpret.motifs.run_tomtom`, which uses memelite's pure-Python `tomtom` (the MEME Suite program is not used). It reads the MEME file, computes p-values, offsets, overlaps and orientation for every query-target pair, and adds an E-value (p-value times the number of targets), a q-value (Benjamini-Hochberg) and consensus strings. Caveat: DNAKit checks that each motif has four columns and its docstring says rows are positions, but gReLU 1.1.0 reads each matrix as shape (4, L) (rows A, C, G, T). The matrix is passed unchanged, so a (length, 4) matrix with length other than 4 is mis-read, and a matrix with four rows and another length is rejected. With four positions, pass the matrix in gReLU's orientation (rows A, C, G, T), as in the example, and check the `Query_consensus` column. tomtom also warns when there are fewer than 25 targets, because its p-values then depend on the target set.
- **API:** `grelu_run_tomtom(motifs: 'Mapping[str, Sequence[Sequence[float]]]', meme_file: 'str', *, verbose: bool = True) -> 'Any'`.
- **Requirements:** Needs the optional packages grelu, memelite and statsmodels (extra `dnakit[bio-grelu]`). It works offline.
- **Input:** `motifs` is a non-empty mapping from names to finite 4-column probability matrices (see the caveat above). `meme_file` is the path of an existing MEME file or one of the bundled names `hocomoco_v12`, `hocomoco_v13` and `consensus`.
- **Output:** A DataFrame with one row per query-target pair and the columns `Query_ID`, `Target_ID`, `Optimal_offset`, `p-value`, `E-value`, `q-value`, `Overlap`, `Query_consensus`, `Target_consensus` and `Orientation`.
- **Example code:**

```python
import warnings

from dnakit.predictions import grelu_run_tomtom

warnings.simplefilter("ignore")
# rows A, C, G, T; four positions: G A T A
query = [
    [0.05, 0.85, 0.05, 0.85],
    [0.05, 0.05, 0.05, 0.05],
    [0.85, 0.05, 0.05, 0.05],
    [0.05, 0.05, 0.85, 0.05],
]
table = grelu_run_tomtom({"query": query}, "db.meme", verbose=False)
best = table.loc[table["p-value"].idxmin()]
print(table.shape, best["Target_ID"], best["Query_consensus"], best["Orientation"])
```

- **Example output:**

```text
(2, 10) GATA_like GATA +
```

## 35) `dnakit.predictions.grelu_variant_to_seqs()`

- **Purpose:** Return the reference and alternate sequences around a single-nucleotide variant.
- **Method:** Calls gReLU `grelu.variant.variant_to_seqs`. It first runs `check_reference` (a mismatch only warns), makes an interval of `seq_len` bases centred on the variant, extracts the reference sequence and inserts the reference and the alternate allele at the centre.
- **API:** `grelu_variant_to_seqs(chrom: 'str', pos: 'int', ref: 'str', alt: 'str', genome: 'str', *, seq_len: 'int' = 1, verbose: bool = True) -> 'tuple[str, str]'`.
- **Requirements:** Needs the optional packages grelu and genomepy (extra `dnakit[bio-grelu]`). A genome name is downloaded on first use; a local FASTA path works offline.
- **Input:** `chrom` and `genome` are non-empty strings; `pos` is the 1-based position; `ref` and `alt` are single bases; `seq_len` is the length of the returned sequences (default 1, which returns only the alleles).
- **Output:** A tuple `(reference_sequence, alternate_sequence)`.
- **Example code:**

```python
from dnakit.predictions import grelu_variant_to_seqs

ref_seq, alt_seq = grelu_variant_to_seqs("chr1", 5, "T", "G", "tiny.fa", seq_len=5, verbose=False)
print(ref_seq, alt_seq)
```

- **Example output:**

```text
GTTGC GTGGC
```

## 36) `dnakit.predictions.grelu_wandb_artifacts()`

- **Purpose:** List the artifacts of a Weights & Biases project of the legacy gReLU model zoo.
- **Method:** Calls `grelu.resources.wandb.artifacts`. It first logs in to Weights & Biases (`wandb.login`, anonymous access) at `host`, then lists the artifact types of the project `grelu/<project>` (optionally those equal to `type_is` or containing `type_contains`) and returns the names of their collections. It fails if no type matches.
- **API:** `grelu_wandb_artifacts(project: 'str', *, host: 'str' = 'https://api.wandb.ai', type_is: 'str | None' = None, type_contains: 'str | None' = None, verbose: bool = True) -> 'list[str]'`.
- **Requirements:** Needs the optional packages grelu and wandb (extra `dnakit[bio-grelu]`) and internet access to Weights & Biases (anonymous access is requested; with an API key you can log in beforehand).
- **Input:** `project` is a non-empty project name of the `grelu` entity; `host` is the W&B API URL; `type_is` and `type_contains` filter the artifact types.
- **Output:** A list of artifact collection names.
- **Example code:**

```python
from dnakit.exceptions import ConfigurationError
from dnakit.predictions import grelu_wandb_artifacts

try:
    grelu_wandb_artifacts("", verbose=False)
except ConfigurationError as error:
    print(type(error).__name__, str(error).split("] ")[1])
```

- **Example output:**

```text
ConfigurationError project must be a non-empty string.
```

Usage with the dependency installed:

```python
from dnakit.predictions import grelu_wandb_artifacts

names = grelu_wandb_artifacts("human-atac-catlas", type_contains="model")
```

## 37) `dnakit.predictions.grelu_wandb_datasets()`

- **Purpose:** List the datasets of a Weights & Biases project of the legacy gReLU model zoo.
- **Method:** Calls `grelu.resources.wandb.datasets`, which is `artifacts` restricted to artifact types containing `dataset`.
- **API:** `grelu_wandb_datasets(project: 'str', *, host: 'str' = 'https://api.wandb.ai', verbose: bool = True) -> 'list[str]'`.
- **Requirements:** Needs the optional packages grelu and wandb (extra `dnakit[bio-grelu]`) and internet access to Weights & Biases.
- **Input:** `project` is a non-empty project name; `host` is the W&B API URL.
- **Output:** A list of dataset artifact names.
- **Example code:**

```python
from dnakit.exceptions import ConfigurationError
from dnakit.predictions import grelu_wandb_datasets

try:
    grelu_wandb_datasets("", verbose=False)
except ConfigurationError as error:
    print(type(error).__name__, str(error).split("] ")[1])
```

- **Example output:**

```text
ConfigurationError project must be a non-empty string.
```

Usage with the dependency installed:

```python
from dnakit.predictions import grelu_wandb_datasets

names = grelu_wandb_datasets("human-atac-catlas")
```

## 38) `dnakit.predictions.grelu_wandb_get_artifact()`

- **Purpose:** Retrieve a Weights & Biases artifact of the legacy gReLU model zoo by name and alias.
- **Method:** Calls `grelu.resources.wandb.get_artifact`, which logs in to Weights & Biases and returns `wandb.Api().artifact("grelu/<project>/<name>:<alias>")`. The artifact object is returned; its files are not downloaded until `artifact.download()` is called.
- **API:** `grelu_wandb_get_artifact(name: 'str', project: 'str', *, host: 'str' = 'https://api.wandb.ai', alias: 'str' = 'latest', verbose: bool = True) -> 'Any'`.
- **Requirements:** Needs the optional packages grelu and wandb (extra `dnakit[bio-grelu]`) and internet access to Weights & Biases.
- **Input:** `name` and `project` are non-empty strings; `alias` is the artifact alias (default `latest`); `host` is the W&B API URL.
- **Output:** A `wandb` artifact object.
- **Example code:**

```python
from dnakit.exceptions import ConfigurationError
from dnakit.predictions import grelu_wandb_get_artifact

try:
    grelu_wandb_get_artifact("", "human-atac-catlas", verbose=False)
except ConfigurationError as error:
    print(type(error).__name__, str(error).split("] ")[1])
```

- **Example output:**

```text
ConfigurationError name must be a non-empty string.
```

Usage with the dependency installed:

```python
from dnakit.predictions import grelu_wandb_get_artifact

artifact = grelu_wandb_get_artifact("model", "human-atac-catlas")
```

## 39) `dnakit.predictions.grelu_wandb_get_dataset_by_model()`

- **Purpose:** List the datasets linked to a model of the legacy gReLU model zoo on Weights & Biases.
- **Method:** Calls `grelu.resources.wandb.get_dataset_by_model`. It gets the model artifact, finds the run that logged it and returns the names of the artifacts that run used.
- **API:** `grelu_wandb_get_dataset_by_model(model_name: 'str', project: 'str', *, host: 'str' = 'https://api.wandb.ai', alias: 'str' = 'latest', verbose: bool = True) -> 'list[str]'`.
- **Requirements:** Needs the optional packages grelu and wandb (extra `dnakit[bio-grelu]`) and internet access to Weights & Biases.
- **Input:** `model_name` and `project` are non-empty strings; `alias` is the model alias; `host` is the W&B API URL.
- **Output:** A list of dataset artifact names.
- **Example code:**

```python
from dnakit.exceptions import ConfigurationError
from dnakit.predictions import grelu_wandb_get_dataset_by_model

try:
    grelu_wandb_get_dataset_by_model("", "human-atac-catlas", verbose=False)
except ConfigurationError as error:
    print(type(error).__name__, str(error).split("] ")[1])
```

- **Example output:**

```text
ConfigurationError model_name must be a non-empty string.
```

Usage with the dependency installed:

```python
from dnakit.predictions import grelu_wandb_get_dataset_by_model

datasets = grelu_wandb_get_dataset_by_model("model", "human-atac-catlas")
```

## 40) `dnakit.predictions.grelu_wandb_get_model_by_dataset()`

- **Purpose:** List the models trained on a dataset of the legacy gReLU model zoo on Weights & Biases.
- **Method:** Calls `grelu.resources.wandb.get_model_by_dataset`. It gets the dataset artifact, takes the first run that used it and returns the names of the artifacts that run logged. It fails if no run used the dataset.
- **API:** `grelu_wandb_get_model_by_dataset(dataset_name: 'str', project: 'str', *, host: 'str' = 'https://api.wandb.ai', alias: 'str' = 'latest', verbose: bool = True) -> 'list[str]'`.
- **Requirements:** Needs the optional packages grelu and wandb (extra `dnakit[bio-grelu]`) and internet access to Weights & Biases.
- **Input:** `dataset_name` and `project` are non-empty strings; `alias` is the dataset alias; `host` is the W&B API URL.
- **Output:** A list of model artifact names.
- **Example code:**

```python
from dnakit.exceptions import ConfigurationError
from dnakit.predictions import grelu_wandb_get_model_by_dataset

try:
    grelu_wandb_get_model_by_dataset("", "human-atac-catlas", verbose=False)
except ConfigurationError as error:
    print(type(error).__name__, str(error).split("] ")[1])
```

- **Example output:**

```text
ConfigurationError dataset_name must be a non-empty string.
```

Usage with the dependency installed:

```python
from dnakit.predictions import grelu_wandb_get_model_by_dataset

models = grelu_wandb_get_model_by_dataset("dataset", "human-atac-catlas")
```

## 41) `dnakit.predictions.grelu_wandb_load_model()`

- **Purpose:** Load a gReLU model from a project of the legacy Weights & Biases model zoo.
- **Method:** Calls `grelu.resources.wandb.load_model`. It gets the model artifact, downloads it to a temporary folder and loads `checkpoint_file` with `LightningModel.load_from_checkpoint`. `weights_only=False` is the default because gReLU checkpoints contain more than tensors; loading a checkpoint from an untrusted source can execute code.
- **API:** `grelu_wandb_load_model(project: 'str', model_name: 'str', *, device: 'str | int' = 'cpu', host: 'str' = 'https://api.wandb.ai', alias: 'str' = 'latest', checkpoint_file: 'str' = 'model.ckpt', weights_only: 'bool' = False, verbose: bool = True) -> 'Any'`.
- **Requirements:** Needs the optional packages grelu and wandb (extra `dnakit[bio-grelu]`), internet access to Weights & Biases and a model checkpoint (downloaded from the artifact).
- **Input:** `project` and `model_name` are non-empty strings; `device` is `cpu` or a GPU index; `alias` is the model alias; `checkpoint_file` is the file name in the artifact; `weights_only` is passed to the checkpoint loader.
- **Output:** A gReLU `LightningModel`.
- **Example code:**

```python
from dnakit.exceptions import ConfigurationError
from dnakit.predictions import grelu_wandb_load_model

try:
    grelu_wandb_load_model("human-atac-catlas", "", verbose=False)
except ConfigurationError as error:
    print(type(error).__name__, str(error).split("] ")[1])
```

- **Example output:**

```text
ConfigurationError model_name must be a non-empty string.
```

Usage with the dependency installed:

```python
from dnakit.predictions import grelu_wandb_load_model

model = grelu_wandb_load_model("human-atac-catlas", "model", device="cpu")
```

## 42) `dnakit.predictions.grelu_wandb_models()`

- **Purpose:** List the models of a Weights & Biases project of the legacy gReLU model zoo.
- **Method:** Calls `grelu.resources.wandb.models`, which is `artifacts` restricted to artifact types containing `model`.
- **API:** `grelu_wandb_models(project: 'str', *, host: 'str' = 'https://api.wandb.ai', verbose: bool = True) -> 'list[str]'`.
- **Requirements:** Needs the optional packages grelu and wandb (extra `dnakit[bio-grelu]`) and internet access to Weights & Biases.
- **Input:** `project` is a non-empty project name; `host` is the W&B API URL.
- **Output:** A list of model artifact names.
- **Example code:**

```python
from dnakit.exceptions import ConfigurationError
from dnakit.predictions import grelu_wandb_models

try:
    grelu_wandb_models("", verbose=False)
except ConfigurationError as error:
    print(type(error).__name__, str(error).split("] ")[1])
```

- **Example output:**

```text
ConfigurationError project must be a non-empty string.
```

Usage with the dependency installed:

```python
from dnakit.predictions import grelu_wandb_models

names = grelu_wandb_models("human-atac-catlas")
```

## 43) `dnakit.predictions.grelu_wandb_projects()`

- **Purpose:** List the projects of the legacy gReLU model zoo on Weights & Biases.
- **Method:** Calls `grelu.resources.wandb.projects`. It logs in to Weights & Biases at `host` (`wandb.login`, anonymous access first) and lists the projects of the `grelu` entity. If the login fails, gReLU raises `RuntimeError`, which DNAKit re-raises as `BackendExecutionError`.
- **API:** `grelu_wandb_projects(*, host: 'str' = 'https://api.wandb.ai', verbose: bool = True) -> 'list[str]'`.
- **Requirements:** Needs the optional packages grelu and wandb (extra `dnakit[bio-grelu]`) and internet access to Weights & Biases.
- **Input:** `host` is the W&B API URL (default `https://api.wandb.ai`).
- **Output:** A list of project names.
- **Example code:**

```python
import os
import tempfile

folder = tempfile.mkdtemp()
os.environ["HOME"] = folder
for name in ("WANDB_CONFIG_DIR", "WANDB_CACHE_DIR", "WANDB_DIR"):
    os.environ[name] = folder
os.environ.pop("WANDB_API_KEY", None)

from dnakit.exceptions import BackendExecutionError
from dnakit.predictions import grelu_wandb_projects

try:
    grelu_wandb_projects(host="http://127.0.0.1:9", verbose=False)
except BackendExecutionError as error:
    print(type(error).__name__, "is not configured" in str(error))
```

- **Example output:**

```text
BackendExecutionError True
```

Usage with the dependency installed:

```python
from dnakit.predictions import grelu_wandb_projects

projects = grelu_wandb_projects()
```

## 44) `dnakit.predictions.grelu_wandb_runs()`

- **Purpose:** List one attribute of every run of a Weights & Biases project of the legacy gReLU model zoo.
- **Method:** Calls `grelu.resources.wandb.runs`, which logs in to Weights & Biases and returns `getattr(run, field)` for every run of `grelu/<project>`, optionally filtered with the `filters` mapping passed to `wandb.Api().runs`.
- **API:** `grelu_wandb_runs(project: 'str', *, host: 'str' = 'https://api.wandb.ai', field: 'str' = 'id', filters: 'Mapping[str, Any] | None' = None, verbose: bool = True) -> 'list[Any]'`.
- **Requirements:** Needs the optional packages grelu and wandb (extra `dnakit[bio-grelu]`) and internet access to Weights & Biases.
- **Input:** `project` is a non-empty project name; `field` is a run attribute such as `id`, `name` or `state`; `filters` is a MongoDB-style filter mapping; `host` is the W&B API URL.
- **Output:** A list with one value per run.
- **Example code:**

```python
from dnakit.exceptions import ConfigurationError
from dnakit.predictions import grelu_wandb_runs

try:
    grelu_wandb_runs("", field="name", verbose=False)
except ConfigurationError as error:
    print(type(error).__name__, str(error).split("] ")[1])
```

- **Example output:**

```text
ConfigurationError project must be a non-empty string.
```

Usage with the dependency installed:

```python
from dnakit.predictions import grelu_wandb_runs

run_ids = grelu_wandb_runs("human-atac-catlas", field="id")
```
