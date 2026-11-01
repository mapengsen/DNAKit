# Kipoi: run, validate and describe local model-zoo models

This group of functions in `dnakit.predictions` is a set of adapters over [kipoi](https://github.com/kipoi/kipoi), the model zoo and API for genomics models. They load kipoi models and dataloaders, describe and list them, run predictions to a TSV file, validate a model against its example input and expected predictions, and run `kipoi test` in the current Python environment or inside a Singularity image. Small helpers expose kipoi's batch handling, example-argument resolution and prediction comparison. The examples and most functions work on local model directories (`source="dir"`); the default source `"kipoi"` clones or updates kipoi's git-backed model zoo and therefore needs internet access.

kipoi is released under the MIT licence (checked in the installed package metadata). The functions are imported with `from dnakit.predictions import ...`. The dependency is optional and is imported only when a function is called; if kipoi is missing, the function raises `BackendUnavailableError` with the install command `python -m pip install kipoi`. The Singularity functions additionally call the external `singularity` program.

All prediction output is written as TSV, because kipoi 0.8.6 cannot write HDF5 files with NumPy 2.

## Example data

All examples use one tiny local model that is also the working directory of the examples: a toy model that returns the GC fraction of one-hot encoded sequences, a dataloader that one-hot encodes the first 20 bases of each FASTA record, and a two-record FASTA file. Where a function needs a zoo (a directory of model folders), the example copies these files into a subfolder of a temporary directory.

**Example file:** `model.yaml`

```text
type: custom
args:
  file: model.py
  object: GcModel
default_dataloader: .
info:
  authors:
    - name: DNAKit docs
  doc: Toy model returning the GC fraction of one-hot encoded sequences.
  version: 0.1.0
  license: MIT
dependencies:
  conda:
    - python
  pip:
    - numpy
schema:
  inputs:
    name: seq
    shape: (None, 4)
    doc: one-hot DNA (A,C,G,T)
  targets:
    name: gc
    shape: (1,)
    doc: GC fraction
```

**Example file:** `model.py`

```text
import numpy as np
from kipoi.model import BaseModel


class GcModel(BaseModel):
    def predict_on_batch(self, x):
        x = np.asarray(x)
        return x[:, :, 1:3].sum(axis=(1, 2)) / x.shape[1]
```

**Example file:** `dataloader.yaml`

```text
type: Dataset
defined_as: dataloader.py::SeqDataset
args:
  fasta_file:
    doc: FASTA file
    example: example.fa
info:
  authors:
    - name: DNAKit docs
  doc: One-hot encodes each FASTA record
output_schema:
  inputs:
    name: seq
    shape: (None, 4)
    doc: one-hot
  metadata:
    name:
      type: str
      doc: record id
dependencies:
  pip:
    - numpy
```

**Example file:** `dataloader.py`

```text
import numpy as np
from kipoi.data import Dataset


class SeqDataset(Dataset):
    def __init__(self, fasta_file):
        self.names, self.seqs = [], []
        for line in open(fasta_file):
            line = line.strip()
            if line.startswith(">"):
                self.names.append(line[1:].split()[0])
                self.seqs.append("")
            elif line:
                self.seqs[-1] += line.upper()

    def __len__(self):
        return len(self.seqs)

    def __getitem__(self, i):
        s = self.seqs[i][:20].ljust(20, "A")
        x = np.zeros((20, 4), dtype="float32")
        for j, c in enumerate(s):
            if c in "ACGT":
                x[j, "ACGT".index(c)] = 1
        return {"inputs": x, "metadata": {"name": self.names[i]}}
```

**Example file:** `example.fa`

```text
>a
GGGGCCCCGGGGCCCCGGGG
>b
AAAAAAAAAATTTTTTTTTT
```

## 1) `dnakit.predictions.compare_kipoi_predictions()`

- **Purpose:** Compare two prediction arrays to a number of decimal places with kipoi's comparison.
- **Method:** Convert both inputs to NumPy arrays and call `kipoi.cli.main.compare_numpy_dict` (the comparison used by `kipoi test`) with `exact=False`. Two arrays match when `abs(expected - observed) < 1.5 * 10**(-decimal)` holds everywhere and the shapes agree. An `AssertionError` from kipoi is turned into `matched=False` and the first line of its message.
- **API:** `compare_kipoi_predictions(observed: 'Any', expected: 'Any', *, decimal: 'int' = 7, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the optional package kipoi (extra `dnakit[bio-kipoi]`).
- **Input:** `observed` and `expected` are numeric arrays or array-likes; nested dictionaries (object arrays) are rejected with `ConfigurationError`. `decimal` is an integer of at least 0.
- **Output:** A dictionary with `matched`, `decimal`, `shape`, `expected_shape`, `mismatch` (the first line of kipoi's message, or `None`) and `provenance`.
- **Example code:**

```python
from dnakit.predictions import compare_kipoi_predictions

same = compare_kipoi_predictions([0.1, 0.2], [0.1, 0.2000001], decimal=5, verbose=False)
print(same["matched"], same["mismatch"])
different = compare_kipoi_predictions([0.1, 0.2], [0.1, 0.3], decimal=5, verbose=False)
print(different["matched"], different["mismatch"])
```

- **Example output:**

```text
True None
False Arrays are not almost equal to 5 decimals
```

## 2) `dnakit.predictions.describe_kipoi_component()`

- **Purpose:** Read the `model.yaml` or `dataloader.yaml` of a local directory and report its description.
- **Method:** Call `kipoi.sources.load_component_descr` on the directory with `which="model"` or `which="dataloader"`. The directory is read as it is and nothing is downloaded.
- **API:** `describe_kipoi_component(component_dir: 'str', *, which: 'str' = 'model', verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the optional package kipoi (extra `dnakit[bio-kipoi]`).
- **Input:** `component_dir` must be an existing directory; `which` is `"model"` or `"dataloader"`, anything else raises `ConfigurationError`.
- **Output:** `type`, `version`, `doc`, `license`, `inputs` and `targets` schemas (each with `name`, `shape` and `doc`), `defined_as` for a dataloader, `component_dir`, `which` and `provenance`. For a model the schemas come from `schema`; for a dataloader from `output_schema`.
- **Example code:**

```python
from dnakit.predictions import describe_kipoi_component

model = describe_kipoi_component(".", verbose=False)
print(model["type"], model["version"], model["license"])
print(model["inputs"]["name"], model["inputs"]["shape"], model["targets"]["name"])
loader = describe_kipoi_component(".", which="dataloader", verbose=False)
print(loader["type"], loader["defined_as"])
```

- **Example output:**

```text
custom 0.1.0 MIT
seq [None, 4] gc
Dataset dataloader.py::SeqDataset
```

## 3) `dnakit.predictions.describe_kipoi_dataloader()`

- **Purpose:** Describe a kipoi dataloader: type, defined class, output schema and arguments.
- **Method:** Call `kipoi.get_dataloader_descr` and report each argument's documentation, type, optional flag and example value. A remote example file is reported as its URL, md5 and name; nothing is downloaded.
- **API:** `describe_kipoi_dataloader(dataloader: 'str', *, source: 'str' = 'kipoi', verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the optional package kipoi (extra `dnakit[bio-kipoi]`); with `source="kipoi"` it also needs internet access to clone or update the model zoo, with `source="dir"` it does not.
- **Input:** `dataloader` is a name in the source or, with `source="dir"`, a directory relative to the current directory. `source` is a kipoi source name or `"dir"`.
- **Output:** `dataloader`, `source`, `provenance` and a `description` with `type`, `defined_as`, `version`, `doc`, `license`, `inputs`, `targets` and an `arguments` mapping.
- **Example code:**

```python
from dnakit.predictions import describe_kipoi_dataloader

result = describe_kipoi_dataloader(".", source="dir", verbose=False)
description = result["description"]
print(description["type"], description["defined_as"])
print(description["inputs"]["shape"])
for name, argument in description["arguments"].items():
    print(name, argument["doc"], argument["optional"], argument["example"])
```

- **Example output:**

```text
Dataset dataloader.py::SeqDataset
[None, 4]
fasta_file FASTA file False example.fa
```

## 4) `dnakit.predictions.get_kipoi_example_kwargs()`

- **Purpose:** Return the example keyword arguments of a model's default dataloader, the information behind `kipoi get-example`.
- **Method:** Load the model description, find its default dataloader description (an imported one is read inside the model directory) and call its `download_example(output_dir=..., absolute_path=True, dry_run=...)`. Remote example files are downloaded into `output_dir` unless `dry_run=True`; plain example values such as `"example.fa"` are returned as written in `dataloader.yaml`.
- **API:** `get_kipoi_example_kwargs(model: 'str', *, source: 'str' = 'kipoi', output_dir: 'str' = 'example', dry_run: 'bool' = False, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the optional package kipoi (extra `dnakit[bio-kipoi]`); the default source and remote example files need internet access, a local `source="dir"` model with plain examples does not.
- **Input:** `model` is a model name or a directory with `source="dir"`; `output_dir` is made absolute; `dry_run` must be a boolean.
- **Output:** `model`, `source`, `output_dir` (absolute), `dry_run`, `example_kwargs` and `provenance`.
- **Example code:**

```python
from pathlib import Path

from dnakit.predictions import get_kipoi_example_kwargs

result = get_kipoi_example_kwargs(".", source="dir", dry_run=True, verbose=False)
print(result["example_kwargs"], result["dry_run"], Path(result["output_dir"]).name)
```

- **Example output:**

```text
{'fasta_file': 'example.fa'} True example
```

## 5) `dnakit.predictions.list_kipoi_directory_dataloaders()`

- **Purpose:** List the dataloaders of a local kipoi model zoo.
- **Method:** Call `kipoi.config.list_dataloaders` on a `LocalSource` built from the directory, so no configured network source is touched. The rows keep the columns `dataloader`, `version`, `type`, `doc`, `inputs`, `targets` and `license`.
- **API:** `list_kipoi_directory_dataloaders(directory: 'str', *, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the optional package kipoi (extra `dnakit[bio-kipoi]`).
- **Input:** `directory` must be an existing directory whose subfolders hold `dataloader.yaml`.
- **Output:** `directory`, `total_dataloaders`, `dataloaders` (a list of row dictionaries) and `provenance`.
- **Example code:**

```python
import shutil
import tempfile
from pathlib import Path

from dnakit.predictions import list_kipoi_directory_dataloaders

zoo = Path(tempfile.mkdtemp())
(zoo / "gc").mkdir()
for name in ["model.yaml", "model.py", "dataloader.yaml", "dataloader.py", "example.fa"]:
    shutil.copy(name, zoo / "gc" / name)
result = list_kipoi_directory_dataloaders(str(zoo), verbose=False)
print(result["total_dataloaders"])
for row in result["dataloaders"]:
    print(row["dataloader"], row["type"], row["version"], row["inputs"])
```

- **Example output:**

```text
1
gc Dataset 0.1 seq
```

## 6) `dnakit.predictions.list_kipoi_directory_models()`

- **Purpose:** List the models of a local kipoi model zoo.
- **Method:** Call `kipoi.config.list_models` on a `LocalSource` built from the directory. The rows keep the columns `model`, `version`, `type`, `doc`, `inputs`, `targets` and `license`.
- **API:** `list_kipoi_directory_models(directory: 'str', *, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the optional package kipoi (extra `dnakit[bio-kipoi]`).
- **Input:** `directory` must be an existing directory whose subfolders hold `model.yaml`; an empty directory gives zero models.
- **Output:** `directory`, `total_models`, `models` (a list of row dictionaries) and `provenance`.
- **Example code:**

```python
import shutil
import tempfile
from pathlib import Path

from dnakit.predictions import list_kipoi_directory_models

zoo = Path(tempfile.mkdtemp())
(zoo / "gc").mkdir()
for name in ["model.yaml", "model.py", "dataloader.yaml", "dataloader.py", "example.fa"]:
    shutil.copy(name, zoo / "gc" / name)
result = list_kipoi_directory_models(str(zoo), verbose=False)
print(result["total_models"])
for row in result["models"]:
    print(row["model"], row["type"], row["version"], row["license"])
```

- **Example output:**

```text
1
gc custom 0.1.0 MIT
```

## 7) `dnakit.predictions.load_kipoi_dataloader()`

- **Purpose:** Return the dataloader class of a kipoi dataloader.
- **Method:** Call `kipoi.get_dataloader_factory`. Instantiate the returned class with the dataloader arguments and iterate it with `batch_iter`.
- **API:** `load_kipoi_dataloader(dataloader: 'str', *, source: 'str' = 'kipoi', verbose: bool = True) -> 'Any'`.
- **Requirements:** Needs the optional package kipoi (extra `dnakit[bio-kipoi]`); the default source needs internet access, `source="dir"` does not.
- **Input:** `dataloader` is a name in the source, or with `source="dir"` a directory relative to the current directory, for example `"my_model/."`. An unknown dataloader raises `BackendExecutionError`.
- **Output:** The dataloader class (a kipoi `Dataset` subclass).
- **Example code:**

```python
from dnakit.predictions import load_kipoi_dataloader

loader_class = load_kipoi_dataloader(".", source="dir", verbose=False)
loader = loader_class(fasta_file="example.fa")
print(loader_class.__name__, len(loader))
for batch in loader.batch_iter(batch_size=2):
    print(sorted(batch), batch["inputs"].shape, batch["metadata"]["name"].tolist())
```

- **Example output:**

```text
SeqDataset 2
['inputs', 'metadata'] (2, 20, 4) ['a', 'b']
```

## 8) `dnakit.predictions.load_kipoi_model()`

- **Purpose:** Load a kipoi model together with its default dataloader.
- **Method:** Call `kipoi.get_model` with the model name, the source and `with_dataloader`. The returned kipoi model object has `predict_on_batch`, `schema`, `info`, `default_dataloader` and the other attributes kipoi attaches.
- **API:** `load_kipoi_model(model: 'str', *, source: 'str' = 'kipoi', with_dataloader: 'bool' = True, verbose: bool = True) -> 'Any'`.
- **Requirements:** Needs the optional package kipoi (extra `dnakit[bio-kipoi]`); the default source clones or updates kipoi's model zoo and needs internet access, `source="dir"` does not download anything.
- **Input:** `model` is a model name or, with `source="dir"`, a directory path (relative to the current directory or absolute). `with_dataloader` must be a boolean.
- **Output:** The kipoi model object.
- **Example code:**

```python
import numpy as np

from dnakit.predictions import load_kipoi_model

model = load_kipoi_model(".", source="dir", verbose=False)
print(type(model).__name__, model.schema.inputs.shape)
x = np.zeros((2, 20, 4), dtype="float32")
x[0, :, 1] = 1
x[1, :, 0] = 1
print(model.predict_on_batch(x).tolist())
```

- **Example output:**

```text
GcModel (None, 4)
[1.0, 0.0]
```

## 9) `dnakit.predictions.predict_kipoi_model()`

- **Purpose:** Run a kipoi model over the data of its default dataloader and write the predictions to a TSV file.
- **Method:** Performs the steps of `kipoi predict`: load the model, instantiate its default dataloader with the validated `dataloader_args`, iterate it with `batch_iter`, call `predict_on_batch` on each batch and write the batches with kipoi's TSV writer after `kipoi.cli.main.prepare_batch`.
- **API:** `predict_kipoi_model(model: 'str', output_file: 'str', *, source: 'str' = 'kipoi', dataloader_args: 'Mapping[str, Any] | None' = None, batch_size: 'int' = 32, num_workers: 'int' = 0, keep_inputs: 'bool' = False, keep_metadata: 'bool' = False, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the optional package kipoi (extra `dnakit[bio-kipoi]`); the default source needs internet access, `source="dir"` does not.
- **Input:** `output_file` must end with `.tsv`, must not exist yet and its folder must exist (HDF5 is refused). `dataloader_args` maps argument names to values; relative paths are relative to the current directory. `batch_size` is at least 1, `num_workers` at least 0; `keep_inputs` and `keep_metadata` are booleans.
- **Output:** `status` (`completed` when the file was written), `model`, `source`, `dataloader` (the class name), `batches`, `samples`, `outputs` (the written path) and `provenance`.
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.predictions import predict_kipoi_model

output = Path(tempfile.mkdtemp()) / "predictions.tsv"
result = predict_kipoi_model(
    ".",
    str(output),
    source="dir",
    dataloader_args={"fasta_file": "example.fa"},
    verbose=False,
)
print(result["status"], result["dataloader"], result["batches"], result["samples"])
print(output.read_text().split())
```

- **Example output:**

```text
completed SeqDataset 1 2
['preds', '1.0', '0.0']
```

## 10) `dnakit.predictions.validate_kipoi_model_example()`

- **Purpose:** Run a local model directory on its example input and optionally compare the predictions with expected ones, as `kipoi test` does.
- **Method:** Load the model with `source="dir"`, call its `pipeline.predict_example` on the dataloader's example arguments, and, when an expectation is available, read the HDF5 file with `kipoi.readers.HDF5Reader` and compare `model.predict_on_batch(inputs)` with the stored `preds` using `compare_numpy_dict` at the precision `test.precision_decimal` of `model.yaml`. The expectation is `expect_file` or a local `test.expect` path in `model.yaml`; a remote `test.expect` is reported in `note` and not downloaded.
- **API:** `validate_kipoi_model_example(model_dir: 'str', *, batch_size: 'int' = 32, output_file: 'str | None' = None, expect_file: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the optional package kipoi (extra `dnakit[bio-kipoi]`; it installs h5py, which reading `expect_file` uses). A model whose example or expectation is remote needs internet access to fetch the example; the local example here does not.
- **Input:** `model_dir` must be an existing directory with `model.yaml`; `output_file` (optional) must be a new `.tsv` path; `expect_file` (optional) must be an existing HDF5 file with `inputs` and `preds` datasets, as written by `kipoi predict -o file.h5 --keep_inputs`, and only a single prediction array is compared.
- **Output:** `status` (`failed` only when the comparison fails), `model_dir`, `samples`, `prediction_shape`, `outputs`, `expect_file`, `compared_samples`, `matched` (`True`, `False` or `None` when nothing was compared), `precision_decimal`, `note` and `provenance`.
- **Example code:**

```python
import tempfile
from pathlib import Path

import h5py
import numpy as np

from dnakit.predictions import validate_kipoi_model_example

plain = validate_kipoi_model_example(".", verbose=False)
print(plain["status"], plain["samples"], plain["prediction_shape"], plain["matched"])

inputs = np.zeros((2, 20, 4), dtype="float32")
inputs[0, :, 1] = 1
inputs[1, :, 0] = 1
expect = Path(tempfile.mkdtemp()) / "expected.h5"
with h5py.File(expect, "w") as handle:
    handle.create_dataset("inputs", data=inputs)
    handle.create_dataset("preds", data=np.array([1.0, 0.0]))
checked = validate_kipoi_model_example(".", expect_file=str(expect), verbose=False)
print(checked["status"], checked["matched"], checked["compared_samples"], checked["precision_decimal"])
with h5py.File(expect, "w") as handle:
    handle.create_dataset("inputs", data=inputs)
    handle.create_dataset("preds", data=np.array([0.0, 1.0]))
wrong = validate_kipoi_model_example(".", expect_file=str(expect), verbose=False)
print(wrong["status"], wrong["matched"])
```

- **Example output:**

```text
completed 2 [2] None
completed True 2 7
failed False
```

## 11) `dnakit.predictions.predict_kipoi_model_in_singularity()`

- **Purpose:** Run `kipoi predict` for a local model directory inside an explicit Singularity image and write TSV predictions.
- **Method:** Build the command `singularity exec -B <dirs> <image> kipoi predict <model_dir> --source dir --batch_size <n> -o <output_file>`, add `--dataloader_args` with the arguments as JSON when they are given, and bind-mount the model directory, the output folder and the folders of local file arguments. The command is run with a timeout of 3600 seconds through the DNAKit external-command runner; a non-zero exit raises `BackendExecutionError` with the log tail. With `dry_run=True` the command is returned and nothing runs.
- **API:** `predict_kipoi_model_in_singularity(container_image: 'str', model_dir: 'str', output_file: 'str', *, dataloader_args: 'Mapping[str, Any] | None' = None, batch_size: 'int' = 32, dry_run: 'bool' = False, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** A real run needs the external program `singularity` (or Apptainer) on PATH and a local `.sif` image that contains kipoi. A missing executable raises `BackendUnavailableError`. The dry run needs neither the program nor kipoi; the image path must exist as a file.
- **Input:** `container_image` must be an existing file, `model_dir` an existing directory, `output_file` a new `.tsv` path in an existing folder; `dataloader_args` is a mapping or `None`; `batch_size` at least 1.
- **Output:** `status` (`dry_run` or `completed`), `operation`, `container_image`, `command`, `bind_directories`, `outputs` (the written TSV on a real run), `return_code`, `log_tail` and `provenance` (the singularity version on a real run).
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.predictions import predict_kipoi_model_in_singularity

folder = Path(tempfile.mkdtemp())
image = folder / "kipoi.sif"
image.write_text("placeholder image file")
result = predict_kipoi_model_in_singularity(
    str(image),
    ".",
    str(folder / "predictions.tsv"),
    dataloader_args={"fasta_file": "example.fa"},
    dry_run=True,
    verbose=False,
)
command = result["command"]
print(result["status"], command[:3], len(command))
print(command[5], command[6], command[-2:])
print(len(result["bind_directories"]), result["outputs"], result["return_code"])
```

- **Example output:**

```text
dry_run ['singularity', 'exec', '-B'] 16
kipoi predict ['--dataloader_args', '{"fasta_file": "example.fa"}']
2 [] None
```

Usage with the dependency installed:

```python
from dnakit.predictions import predict_kipoi_model_in_singularity

result = predict_kipoi_model_in_singularity(
    "kipoi.sif",
    "my_model",
    "predictions.tsv",
    dataloader_args={"fasta_file": "reads.fa"},
)
print(result["status"], result["outputs"])
```

## 12) `dnakit.predictions.test_kipoi_model_in_singularity()`

- **Purpose:** Run `kipoi test` on a local model directory inside an explicit Singularity image.
- **Method:** Build and run `singularity exec -B <model dir> <image> kipoi test --batch_size <n> --source dir <model_dir>`, the same `kipoi test` that `kipoi test-source --singularity` runs for each model. kipoi's own test compares the predictions with a local `test.expect`. The run has a timeout of 1800 seconds. A non-zero exit raises `BackendExecutionError` with the log tail; with `dry_run=True` the command is returned instead.
- **API:** `test_kipoi_model_in_singularity(container_image: 'str', model_dir: 'str', *, batch_size: 'int' = 32, dry_run: 'bool' = False, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** A real run needs the external program `singularity` (or Apptainer) on PATH and a local `.sif` image that contains kipoi. This environment has no singularity, so the example below shows the error raised when the program is missing.
- **Input:** `container_image` must be an existing file, `model_dir` an existing directory, `batch_size` at least 1, `dry_run` a boolean.
- **Output:** The same dictionary as the predict function: `status`, `operation`, `container_image`, `command`, `bind_directories`, `outputs` (empty), `return_code`, `log_tail` and `provenance`.
- **Example code:**

```python
import os
import tempfile
from pathlib import Path

from dnakit.exceptions import BackendUnavailableError
from dnakit.predictions import test_kipoi_model_in_singularity

image = Path(tempfile.mkdtemp()) / "kipoi.sif"
image.write_text("placeholder image file")
os.environ["PATH"] = ""
try:
    test_kipoi_model_in_singularity(str(image), ".", verbose=False)
except BackendUnavailableError as error:
    print(type(error).__name__, "singularity executable was not found" in str(error))
```

- **Example output:**

```text
BackendUnavailableError True
```

Usage with the dependency installed:

```python
from dnakit.predictions import test_kipoi_model_in_singularity

result = test_kipoi_model_in_singularity("kipoi.sif", "my_model", batch_size=4)
print(result["status"], result["return_code"])
```

## 13) `dnakit.predictions.test_kipoi_model_in_current_env()`

- **Purpose:** Run `kipoi test` on a local model directory in the current Python environment and fail on logged warnings.
- **Method:** Call kipoi's `cli_test` with `--source dir --batch_size <n>` on the directory, with stdout and stderr captured. The example inputs from `dataloader.yaml` are predicted and the local `test.expect` of `model.yaml`, if any, is compared at the model's precision. Any `WARNING` record logged during the test makes the call raise `BackendExecutionError`, as `kipoi.cli.source_test.test_model` does for its subprocess; a kipoi exception or non-zero exit also raises it. Nothing is installed and no conda environment is created.
- **API:** `test_kipoi_model_in_current_env(model_dir: 'str', *, batch_size: 'int' = 32, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the optional package kipoi (extra `dnakit[bio-kipoi]`) and the model's own dependencies already importable.
- **Input:** `model_dir` must be an existing directory; `batch_size` at least 1.
- **Output:** `status`, `operation`, `model_dir`, `batch_size`, `warnings` (always 0 on success), `log_tail` (the captured output, which can include a progress bar) and `provenance`.
- **Example code:**

```python
from pathlib import Path

from dnakit.predictions import test_kipoi_model_in_current_env

result = test_kipoi_model_in_current_env(".", batch_size=2, verbose=False)
print(result["status"], result["batch_size"], result["warnings"], Path(result["model_dir"]).is_dir())
```

- **Example output:**

```text
completed 2 0 True
```

## 14) `dnakit.predictions.prepare_kipoi_batch()`

- **Purpose:** Attach model predictions to one dataloader batch the way `kipoi predict` does.
- **Method:** Copy the batch into a new dictionary and call `kipoi.cli.main.prepare_batch`. The predictions are stored under `"preds"`. Unless `keep_inputs` is true, the `"inputs"` and `"targets"` entries are removed; unless `keep_metadata` is true, `"metadata"` is removed. The caller's dictionary is not changed.
- **API:** `prepare_kipoi_batch(dataloader_batch: 'dict[str, Any]', predictions: 'Any', *, keep_inputs: 'bool' = False, keep_metadata: 'bool' = False, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the optional package kipoi (extra `dnakit[bio-kipoi]`).
- **Input:** `dataloader_batch` is a mapping such as one batch yielded by a kipoi dataloader (a non-mapping raises `ConfigurationError`); `predictions` is the model output; `keep_inputs` and `keep_metadata` are booleans.
- **Output:** A new dictionary with `preds` and, depending on the flags, `inputs`, `targets` and `metadata`.
- **Example code:**

```python
import numpy as np

from dnakit.predictions import prepare_kipoi_batch

batch = {
    "inputs": np.zeros((2, 3)),
    "targets": np.ones(2),
    "metadata": {"name": np.array(["a", "b"])},
}
predictions = np.array([0.5, 0.25])
print(sorted(prepare_kipoi_batch(batch, predictions, verbose=False)))
full = prepare_kipoi_batch(batch, predictions, keep_inputs=True, keep_metadata=True, verbose=False)
print(sorted(full), full["preds"].tolist())
print(sorted(batch))
```

- **Example output:**

```text
['preds']
['inputs', 'metadata', 'preds', 'targets'] [0.5, 0.25]
['inputs', 'metadata', 'targets']
```

## 15) `dnakit.predictions.concatenate_kipoi_batches()`

- **Purpose:** Join dataloader batches along the first axis, field by field.
- **Method:** Call `kipoi.cli.main.numpy_collate_concat` on the list of batches: arrays are joined with `numpy.concatenate`, and nested dictionaries or lists are joined entry by entry. The result has the structure of one batch with all the samples of all batches.
- **API:** `concatenate_kipoi_batches(batches: 'Sequence[Any]', *, verbose: bool = True) -> 'Any'`.
- **Requirements:** Needs the optional package kipoi (extra `dnakit[bio-kipoi]`).
- **Input:** A non-empty sequence of batches with the same structure. A string, bytes object or mapping, or an empty sequence, raises `ConfigurationError`.
- **Output:** One batch with the same structure, holding the samples of all input batches in order.
- **Example code:**

```python
import numpy as np

from dnakit.predictions import concatenate_kipoi_batches

first = {"inputs": np.zeros((2, 3)), "metadata": {"name": np.array(["a", "b"])}}
second = {"inputs": np.ones((3, 3)), "metadata": {"name": np.array(["c", "d", "e"])}}
joined = concatenate_kipoi_batches([first, second], verbose=False)
print(joined["inputs"].shape, joined["inputs"][:, 0].tolist())
print(joined["metadata"]["name"].tolist())
```

- **Example output:**

```text
(5, 3) [0.0, 0.0, 1.0, 1.0, 1.0]
['a', 'b', 'c', 'd', 'e']
```

## 16) `dnakit.predictions.kipoi_dataloader_example_kwargs()`

- **Purpose:** Resolve the example arguments of a local kipoi dataloader.
- **Method:** Load `dataloader.yaml` with `kipoi.sources.load_component_descr` and call `kipoi.specs.example_kwargs` on its arguments. Plain values are returned as written; a remote example file becomes its local path under `cache_dir` (downloaded there when missing or when its md5 differs, unless `dry_run`), or, without `cache_dir`, a dictionary with its `url`, `md5` and `name`. Nothing is downloaded in the examples because `dry_run` or no `cache_dir` is used.
- **API:** `kipoi_dataloader_example_kwargs(dataloader_dir: 'str', *, cache_dir: 'str | None' = None, absolute_path: 'bool' = True, dry_run: 'bool' = False, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the optional package kipoi (extra `dnakit[bio-kipoi]`); downloading remote example files (a `cache_dir` without `dry_run`) needs internet access.
- **Input:** `dataloader_dir` must be an existing directory holding `dataloader.yaml`; `cache_dir` is optional; `absolute_path` and `dry_run` are booleans.
- **Output:** `dataloader_dir`, `cache_dir`, `dry_run`, `example_kwargs` and `provenance`. With `cache_dir`, the cache folder is created even when `dry_run` is true.
- **Example code:**

```python
import tempfile
import textwrap
from pathlib import Path

from dnakit.predictions import kipoi_dataloader_example_kwargs

print(kipoi_dataloader_example_kwargs(".", verbose=False)["example_kwargs"])

folder = Path(tempfile.mkdtemp())
(folder / "dataloader.py").write_text("from kipoi.data import Dataset\n\n\nclass Loader(Dataset):\n    pass\n")
(folder / "dataloader.yaml").write_text(textwrap.dedent("""\
    type: Dataset
    defined_as: dataloader.py::Loader
    args:
      fasta_file:
        doc: FASTA file
        example:
          url: https://example.org/data/example.fa
          md5: d41d8cd98f00b204e9800998ecf8427e
      window:
        doc: window size
        example: 10
    info:
      authors:
        - name: DNAKit docs
      doc: Remote example
    output_schema:
      inputs:
        name: seq
        shape: (None, 4)
        doc: one-hot
    """))
remote = kipoi_dataloader_example_kwargs(str(folder), verbose=False)["example_kwargs"]
print(remote["fasta_file"]["url"], remote["window"])
cache = folder / "cache"
planned = kipoi_dataloader_example_kwargs(str(folder), cache_dir=str(cache), dry_run=True, verbose=False)
print(Path(planned["example_kwargs"]["fasta_file"]).relative_to(cache), cache.exists())
```

- **Example output:**

```text
{'fasta_file': 'example.fa'}
https://example.org/data/example.fa 10
fasta_file True
```

## 17) `dnakit.predictions.kipoi_dataloader_default_args()`

- **Purpose:** Resolve the default arguments of a local kipoi dataloader.
- **Method:** Load `dataloader.yaml` and call `kipoi.specs.download_default_args`, the step kipoi runs when it loads a model. Remote default files are downloaded into `output_dir` and replaced by their local paths; other non-`None` defaults are returned as written, existing file paths are made absolute and booleans stay flags. Numeric, list, tuple and dictionary defaults bypass the upstream call, because kipoi 0.8.6 passes every non-boolean default to `os.path.exists`, which reads an integer as a file descriptor.
- **API:** `kipoi_dataloader_default_args(dataloader_dir: 'str', output_dir: 'str', *, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the optional package kipoi (extra `dnakit[bio-kipoi]`); remote default files need internet access, a dataloader without remote defaults downloads nothing.
- **Input:** `dataloader_dir` must be an existing directory holding `dataloader.yaml`; `output_dir` is made absolute.
- **Output:** `dataloader_dir`, `output_dir`, `overrides` (a mapping of argument name to resolved default) and `provenance`. Arguments with no default are absent. Note that the numeric, list and dictionary defaults are listed after the other defaults.
- **Example code:**

```python
import tempfile
import textwrap
from pathlib import Path

from dnakit.predictions import kipoi_dataloader_default_args

folder = Path(tempfile.mkdtemp())
(folder / "dataloader.py").write_text("from kipoi.data import Dataset\n\n\nclass Loader(Dataset):\n    pass\n")
(folder / "dataloader.yaml").write_text(textwrap.dedent("""\
    type: Dataset
    defined_as: dataloader.py::Loader
    args:
      fasta_file:
        doc: FASTA file
      window:
        doc: window size
        default: 20
      strict:
        doc: switch
        default: true
      label:
        doc: label
        default: x
    info:
      authors:
        - name: DNAKit docs
      doc: Defaults
    output_schema:
      inputs:
        name: seq
        shape: (None, 4)
        doc: one-hot
    """))
output = folder / "downloads"
result = kipoi_dataloader_default_args(str(folder), str(output), verbose=False)
print(dict(sorted(result["overrides"].items())), output.exists())
```

- **Example output:**

```text
{'label': 'x', 'strict': True, 'window': 20} False
```

## 18) `dnakit.predictions.test_kipoi_source_models()`

- **Purpose:** Run `kipoi test` on the models of a local model source in the current Python environment.
- **Method:** List the model folders of `source_dir` (or validate the names in `models`), then run `test_kipoi_model_in_current_env` for each, so a logged `WARNING` counts as a failure. Failures are collected in the result instead of stopping the run, unless `stop_on_failure` is true. This is the prediction part of `kipoi test-source`; conda environments and git ranges are not used.
- **API:** `test_kipoi_source_models(source_dir: 'str', *, models: 'Sequence[str] | None' = None, batch_size: 'int' = 4, stop_on_failure: 'bool' = False, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the optional package kipoi (extra `dnakit[bio-kipoi]`) and the models' own dependencies already importable.
- **Input:** `source_dir` must be an existing directory; `models` is an optional sequence of subfolder names (absolute names, `..` and non-folders raise `ConfigurationError`, as does a directory with no models); `batch_size` at least 1; `stop_on_failure` a boolean.
- **Output:** `status` (`passed` or `failed`), `source_dir`, `batch_size`, `tested`, `passed`, `failed` (the failed names), `results` (one row per model with `model`, `status`, `warnings` and `error`) and `provenance`.
- **Example code:**

```python
import shutil
import tempfile
from pathlib import Path

from dnakit.predictions import test_kipoi_source_models

zoo = Path(tempfile.mkdtemp())
for folder in ["gc_a", "gc_b"]:
    (zoo / folder).mkdir()
    for name in ["model.yaml", "model.py", "dataloader.yaml", "dataloader.py", "example.fa"]:
        shutil.copy(name, zoo / folder / name)
result = test_kipoi_source_models(str(zoo), verbose=False)
print(result["status"], result["tested"], result["passed"], result["failed"])
for row in result["results"]:
    print(row["model"], row["status"], row["warnings"], row["error"])
```

- **Example output:**

```text
passed 2 2 []
gc_a passed 0 None
gc_b passed 0 None
```
