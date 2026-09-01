# Engineering and scaling

Unified management of DNAKit's computing backend, Python API, CLI, configuration workflow, batch and parallel computing, cache, version and recurrence information.

## 1) `ENG-003` Python API

- **What it does:** Call DNAKit functions through stable Python objects, functions, and result types, making it easy to combine analysis workflows in scripts, notebooks, tests, and other software.
- **API**: The `dnakit` top-level and `dnakit.<domain>` modules are namespace entries and have no direct calling parameters; see the corresponding function entries for each function parameter.
- **Input**: `DNASequence`, `DNARecord`, `DNASet` or configuration object required by the corresponding function.
- **Sample Code**:

```python
import dnakit
from dnakit.descriptors import length_features

sequence = dnakit.normalize(" acgt ").sequence
assert sequence is not None
print(length_features(sequence).symbol_length)  # 4
```

- **Example results:**

```text
4
```

## 2) `ENG-004` CLI

- **Function:** Run common functions such as standardization, description, fingerprinting, search, and comparison through command line parameters, and output structured JSON, suitable for shell processes and automation tasks.
- **API**: `dnakit COMMAND[required] ARGS[optional]`; the complete parameters of each subcommand are subject to `dnakit COMMAND --help`.
- **Input**: subcommand and its sequence, file or configuration parameters.
- **Sample Code**:

```bash
dnakit describe ACGTACGT
dnakit compare ACGT ACGA --method hamming
```

- **Example results:**

```text
{"base_composition": {"ambiguity_policy": "ignore", "counts": {"A": 2, "C": 2, "G": 2, "T": 2}, "cross_gaps": false, "denominator": 8, "fractions": {"A": 0.25, "C": 0.25, "G": 0.25, "T": 0.25}, "gap_count": 0, "ignored_ambiguity_count": 0, "method": "canonical_base_count", "name": "base_composition", "sequence_id": null, "unknown_gap_count": 0}, "complexity": {"ambiguity_policy": "ignore", "by_k": {"1": 1.0, "2": 0.5714285714285714, "3": 0.6666666666666666, "4": 0.8, "5": 1.0, "6": 1.0}, "cross_gaps": false, "formula": "product_k(unique_kmers/min(4**k,valid_kmer_positions))", "gap_count": 0, "max_observations": 10000000, "max_word_size": 6, "method": "vocabulary-observed-over-possible-product", "name": "linguistic_complexity", "observation_count": 33, "observed_by_k": {"1": 4, "2": 4, "3": 4, "4": 4, "5": 4, "6": 3}, "possible_by_k": {"1": 4, "2": 7, "3": 6, "4": 5, "5": 4, "6": 3}, "score": 0.3047619047619048, "sequence_id": null, "unknown_gap_count": 0}, "gc_at": {"ambiguity_policy": "ignore", "at_count": 4, "at_fraction": 0.5, "cross_gaps": false, "denominator": 8, "gap_count": 0, "gc_count": 4, "gc_fraction": 0.5, "ignored_ambiguity_count": 0, "method": "canonical_base_fraction", "name": "gc_at_content", "sequence_id": null, "unknown_gap_count": 0}, "repeat": {"ambiguity_policy": "ignore", "comparisons": 4, "cross_gaps": false, "denominator": 8, "gap_count": 0, "max_comparisons": 5000000, "max_unit_length": 20, "method": "maximal-exact-tandem-repeat-union", "min_repeats": 2, "min_unit_length": 1, "name": "exact_repeat_fraction", "repeat_count_by_unit": {"4": 1}, "repeat_fraction": 1.0, "repeated_base_count": 8, "runs": [{"repeat_count": 2, "symbol_end": 8, "symbol_start": 0, "unit": "ACGT", "unit_length": 4}], "sequence_id": null, "unknown_gap_count": 0}}
{"costs": {"substitution": 1.0}, "distance": 1.0, "dp_cells": null, "edit_path": null, "exceeded_max_distance": false, "iupac_matching": "literal", "left_id": null, "left_length": 4, "max_cells": null, "max_distance": null, "method": "hamming", "mismatches": [{"left_symbol": "T", "position": 3, "right_symbol": "A"}], "name": "hamming_distance", "right_id": null, "right_length": 4}
```

## 3) `ENG-006` Batch calculation

- **Use:** Repeatedly call the same calculation on a large number of DNA records, collecting successful results and errors one by one in the order of input, avoiding a single failure interrupting the entire batch.
- **API**: `dnakit.batch.run_batch(records[required], operation[required], name[required], config[optional], progress[optional])`, `dnakit.batch.iter_batch(records[required], operation[required], config[optional], progress[optional])`; `config` uses `dnakit.batch.BatchConfig`.
- **Input**: Record iterator, batch callable and optional error, seed, resume configuration.
- **Sample Code**:

```python
from dnakit import DNARecord, DNASequence
from dnakit.batch import run_batch

records = [
    DNARecord(DNASequence("A"), "a"),
    DNARecord(DNASequence("CC"), "b"),
]
result = run_batch(
    records,
    lambda record, context: record.sequence.symbol_length,
    name="length",
)
print([item.value for item in result.items])  # [1, 2]
```

- **Example results:**

```text
[1, 2]
```

## 4) `ENG-007` Parallel calculation

- **Function:** Use controlled worker threads to execute independent recording tasks in parallel, while keeping the results consistent with the input order, and uniformly propagating cancellation, progress and error information.
- **API**: `dnakit.batch.run_batch(records[required], operation[required], name[required], config[optional], progress[optional])`; `config` uses `dnakit.batch.BatchConfig`, this setting is `execution_mode="thread"`.
- **Input**: record, callable, number of workers and optional `max_in_flight`.
- **Sample Code**:

```python
from dnakit import DNARecord, DNASequence
from dnakit.batch import BatchConfig, run_batch

records = [DNARecord(DNASequence("ACGT"), str(index)) for index in range(4)]
result = run_batch(
    records,
    lambda record, context: (record.id, context.seed),
    name="threaded",
    config=BatchConfig(seed=7, jobs=2, execution_mode="thread", max_in_flight=2),
)
print(result.success_count)  # 4
```

- **Example results:**

```text
4
```

## 5) `ENG-008` Chunking and streaming

- **Function:** Read, convert and write large files in blocks using iterators, so that only the current batch is retained in the memory, and progress events are returned, which is suitable for data that exceeds the memory capacity.
- **API**: `dnakit.read(source[required], format[optional], config[optional])`, `dnakit.RecordSource(iterator[required], close_callback[optional], source_name[optional], format[optional])`, `dnakit.io.iter_chunks(values[required], chunk_size[required])`, `dnakit convert input_path[required] output_path[required] --input-format[optional] --output-format[optional] --overwrite[optional] --progress/--no-progress[optional]`.
- **Input**: file path or iterator; optional format, compression and chunk size.
- **Sample Code**:

```python
from dnakit.io import iter_chunks

for chunk in iter_chunks(range(5), chunk_size=2):
    print(chunk)
# (0, 1), (2, 3), (4,)
```

- **Example results:**

```text
(0, 1)
(2, 3)
(4,)
```

## 6) `ENG-009` cache

- **Function:** Generate cache keys based on input content, function names and parameters, reuse existing calculation results and verify integrity, reducing repeated high-cost calculations.
- **API**: `dnakit.cache.CacheKey(namespace[required], digest[required], schema_version[optional])`, `dnakit.cache.CacheKey.from_components(namespace[required], components[required], schema_version[optional])`, `dnakit.cache.JSONCache(root[required], max_entry_bytes[optional])`, `dnakit.cache.JSONCache.get(key[required])`, `dnakit.cache.JSONCache.put(key[required], value[required])`, `dnakit.cache.JSONCache.invalidate(key[required])`, `dnakit.cache.JSONCache.clear(namespace[optional])`.
- **Input**: Private cache directory, namespace, input/parameters/algorithm and backend version components.
- **Sample Code**:

```python
from pathlib import Path
from tempfile import TemporaryDirectory

from dnakit.cache import CacheKey, JSONCache

with TemporaryDirectory() as directory:
    cache = JSONCache(Path(directory) / "cache")
    key = CacheKey.from_components("length", {"sequence": "ACGT", "algorithm": "v1"})
    cache.put(key, {"value": 4})
    print(cache.get(key))  # {"value": 4}
```

- **Example results:**

```text
{'value': 4}
```

## 7) `ENG-012` Errors and Warnings

- **Function:** Use structured exceptions and Issues to uniformly return error codes, context, severity and repair tips, so that CLI, API and reports can consistently explain failures or data problems.
- **API**: `dnakit.exceptions.DNAKitError(message[required], code[optional], context[optional], hint[optional])`, `dnakit.core.Issue(code[required], severity[required], message[required], location[optional], details[optional])`; CLI non-zero exit code without calling parameters.
- **Input**: Illegal input, misconfigured or missing backend request.
- **Sample Code**:

```python
from dnakit import DNASequence
from dnakit.exceptions import DNAKitError

try:
    DNASequence("AX")
except DNAKitError as error:
    print(error.code)     # INVALID_ALPHABET
    print(error.context)
```

- **Example results:**

```text
INVALID_ALPHABET
{'alphabet': 'strict', 'part_index': 0, 'part_offset': 1, 'symbol': 'X'}
```

## 8) `ENG-013` Agents and MCP

- **Function:** Convert stable public functions into searchable, inspectable, and callable MCP tools while reusing the existing DNAKit implementations.
- **API**: `dnakit-mcp`, `dnakit.tools.default_tool_registry()`, and `dnakit.tools.create_server()`.
- **Input**: A tool name and JSON arguments matching the generated schema; file writes, model downloads, and external programs also require explicit authorization.
- **Sample Code**:

```python
from dnakit.tools import default_tool_registry

registry = default_tool_registry()
result = registry.execute(
    "dnakit.thermodynamics.molecular_weight",
    {"sequence": "ACGT"},
)
print(result["value_dalton"])
```

- **Example result:**

```text
1173.84
```

See [Agent and MCP tools](../../agent_tools.md) for client configuration.
