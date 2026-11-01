# P41 sequence operations

## Purpose

Complete selected P41 operations through dedicated upstream APIs, returning JSON data or actual file paths.

## Method

Upstream algorithms and defaults are preserved. Choose algorithm/output variants with explicit parameters. Existing outputs are refused. Optional dependencies are required only when called.

## API

| Module | Functions |
| --- | --- |
| `dnakit.alignment.edlib_p41_functions` | `edlib_prefix_alignment` |
| `dnakit.annotation.ncls_p41_functions` | `ncls_limited_overlaps` |
| `dnakit.phylogeny.pyvolve_p41_functions` | `simulate_pyvolve_sequences`, `pyvolve_state_frequencies` |
| `dnakit.regulation.tangermeme_p41_functions` | `tangermeme_hypothetical_attributions`, `tangermeme_unchunk`, `tangermeme_extract_signal`, `tangermeme_position_entropy` |
| `dnakit.sequence.py2bit_p41_functions` | `twobit_region_composition`, `twobit_masked_blocks` |
| `dnakit.sequence.screed_p41_functions` | `screed_record_at_index`, `screed_record_names`, `screed_is_dna`, `screed_write_record` |
| `dnakit.sequencing_qc.pyslow5_p41_functions` | `slow5_read_records` |

## Input

Use explicit local file paths, JSON-compatible arrays and 0-based half-open intervals unless an API documents a native alternative. External CLI programs must be installed. Refer to each function docstring for limits and parameter choices.

## Output

Dictionaries contain status and parsed results, or paths of files actually written. Plot interfaces also return the input sequence and features. Tests and limitations: `paper/write_paper/figure/fig2/task_coverage/coverage_gap/manifests/integration_p41.json`.

### SHW alignment

**Example code:**

```python
from dnakit.alignment.edlib_p41_functions import edlib_prefix_alignment
r = edlib_prefix_alignment('ACGT', 'ACGTAAA', task='path', verbose=False)
print(r['edit_distance'], r['cigar'])
```

**Example output:**

```text
0 4=
```

### Limited interval overlaps

**Example code:**

```python
from dnakit.annotation.ncls_p41_functions import ncls_limited_overlaps
r = ncls_limited_overlaps([0, 5], [4, 9], [1], [3], selection='first', verbose=False)
print(r['query_indices'], r['stored_indices'])
```

**Example output:**

```text
[0] [0]
```

### FASTX write and alphabet

**Example code:**

```python
from tempfile import TemporaryDirectory
from pathlib import Path
from dnakit.sequence.screed_p41_functions import screed_write_record, screed_is_dna
with TemporaryDirectory() as folder:
    p = Path(folder) / 'read.fa'
    screed_write_record({'name': 'read', 'sequence': 'ACGTN'}, str(p), verbose=False)
    print(p.read_text().strip())
print(screed_is_dna('ACGTN', verbose=False)['is_dna'])
```

**Example output:**

```text
>read
ACGTN
True
```

### Nucleotide frequencies

**Example code:**

```python
from tempfile import TemporaryDirectory
from pathlib import Path
from dnakit.phylogeny.pyvolve_p41_functions import pyvolve_state_frequencies
with TemporaryDirectory() as folder:
    p = Path(folder) / 'dna.fa'
    p.write_text('>dna\nACGT\n')
    print(pyvolve_state_frequencies(str(p), verbose=False)['frequencies'])
```

**Example output:**

```text
{'A': 0.25, 'C': 0.25, 'G': 0.25, 'T': 0.25}
```

### Entropy and chunk reconstruction

**Example code:**

```python
from dnakit.regulation.tangermeme_p41_functions import tangermeme_position_entropy, tangermeme_unchunk
r = tangermeme_position_entropy([[[0.25], [0.25], [0.25], [0.25]]], verbose=False)
print(round(r['entropy'][0][0], 6))
print(tangermeme_unchunk([[[1, 2]], [[3, 4]]], [4], verbose=False)['sequences'])
```

**Example output:**

```text
2.0
[[[1.0, 2.0, 3.0, 4.0]]]
```

### SLOW5 read IDs

**Example code:**

```python
from tempfile import TemporaryDirectory
from pathlib import Path
from dnakit.sequencing_qc.backends import write_slow5_signals
from dnakit.sequencing_qc.pyslow5_p41_functions import slow5_read_records
with TemporaryDirectory() as folder:
    p = Path(folder) / 'read.slow5'
    write_slow5_signals([{'read_id':'r1','signal':[1,2,3],'digitisation':8192,'offset':0,'range':1467,'sampling_rate':4000}], p, verbose=False)
    print(slow5_read_records(str(p), include_signal=False, verbose=False)['records'])
```

**Example output:**

```text
['r1']
```

### 2bit masked blocks

**Example code:**

```python
from pathlib import Path
import py2bit
from dnakit.sequence.py2bit_p41_functions import twobit_masked_blocks
# This small fixture is distributed with py2bit; supply your own .2bit in production.
p = Path(py2bit.__file__).parent / 'py2bitTest' / 'foo.2bit'
print(twobit_masked_blocks(str(p), 'chr1', kind='hard', verbose=False)['blocks'])
```

**Example output:**

```text
[[0, 50], [100, 150]]
```

