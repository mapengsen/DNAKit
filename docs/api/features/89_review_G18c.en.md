# G18c review: sbol3, screed and pyfaidx

This page records the G18c review. The ledger has 123 rows with status `other`: 118 in sbol3, 3 in screed and 2 in pyfaidx. After the review every row has a final status and a one-line reason. This group adds no new wrapper. The three modules `sbol3_review_functions.py`, `screed_review_functions.py` and `pyfaidx_review_functions.py` export no function; they are kept for later wiring.

- No new DNA-level function in this group runs offline on local input. The 118 sbol3 rows are global namespace configuration helpers, class-module re-exports of them, and a package-path helper. The 5 screed and pyfaidx rows are the HAVA example format, two file-name and modification-time utilities, and one pair writer that duplicates an existing wrapper.
- The existing DNA-level wrappers stay where they are: sbol3 `copy_sbol3_design` and `to_sbol3_display_id` in `molecular_lab/sbol3_functions.py`; screed and pyfaidx in `sequence/screed_functions.py` and `sequence/pyfaidx_functions.py`.
- The test file `others/tests/unit/test_review_G18c.py` calls real upstream APIs or works hand calculations to check each classification.

## Review results

| Package | Other rows before | New status | New wrapped |
| --- | --- | --- | --- |
| sbol3 | 118 | internal 118 | 0 |
| screed | 3 | not-dna 2, duplicate-of 1 (`format_fastx_records`) | 0 |
| pyfaidx | 2 | not-dna 2 | 0 |

- sbol3, 118 rows: every class module re-exports pySBOL3's `get_namespace`, `set_namespace` and `set_defaults`, which is 38 modules × 3 = 114 rows. The top-level functions with the same names add 3 rows, and `sbol3.document.data_path` adds 1. The namespace functions read and write the process-wide namespace, and `data_path` only joins a relative path onto the sbol3 install directory. Under the strict rule all 118 are marked `internal` and not wrapped.
- screed `hava_iter` and `read_hava_sequences`: upstream HAVA is an example format for the API reference. Each record is six fixed lines with placeholder field names (`hava`, `quarzk` and so on) and has no DNA meaning, so both are `not-dna`.
- screed `write_fastx_pair`: upstream calls the undefined name `write_record` and raises NameError. Its intended output is the same two records written in order, which `format_fastx_records` already produces, so it is `duplicate-of:format_fastx_records`.
- pyfaidx `getmtime_fsspec` reads file modification times and `get_valid_filename` sanitizes file names. Neither touches sequence data, so both are `not-dna`.

## 1) `format_fastx_records` covers `write_fastx_pair`

- **Purpose:** Write two records in order as FASTA or FASTQ text. This covers the purpose of upstream `screed.screedRecord.write_fastx_pair`, which cannot run upstream.
- **Method:** `format_fastx_records` calls `screed.screedRecord.write_fastx` once per record, so its output equals two single-record writes. In FASTQ mode every record must carry a quality string of the same length as its sequence; upstream only checks that both records carry a quality.
- **API:** `format_fastx_records(records[required], output_format[optional], max_records[optional], verbose[optional])`
- **Input:** A list of mappings, each with `name` and `sequence`; FASTQ also needs `quality`. `description` is optional.
- **Output:** A dictionary with `record_count`, `text` and `provenance`.

**Example code:**

```python
from dnakit.sequence import format_fastx_records

pair = [
    {"name": "read1", "sequence": "ACGT"},
    {"name": "read2", "sequence": "GGCC"},
]
fasta = format_fastx_records(pair, output_format="fasta", verbose=False)
print(fasta["record_count"])
print(fasta["text"], end="")

fastq = format_fastx_records(
    [dict(pair[0], quality="IIII")],
    output_format="fastq",
    verbose=False,
)
print(fastq["text"], end="")
```

**Example output:**

```text
2
>read1
ACGT
>read2
GGCC
@read1
ACGT
+
IIII
```

## 2) Why the sbol3 namespace helpers are not wrapped

- **Purpose:** Explain why `sbol3.component.set_namespace`, `sbol3.attachment.set_defaults` and the other 114 class-module re-exports are not wrapped. They are the same objects as the top-level `sbol3.set_namespace` and `sbol3.set_defaults`.
- **Method:** These functions read and write pySBOL3's process-wide namespace. `set_namespace` accepts only a URL and raises ValueError otherwise; `set_defaults` resets the namespace to None. A wrapper would change every later SBOL3 call in the process, so none is provided.
- **API:** None; this group exports nothing.
- **Input:** None.
- **Output:** None.

**Example code:**

```python
import sbol3

sbol3.component.set_namespace("https://example.org/lab")
print(sbol3.get_namespace())
sbol3.attachment.set_defaults()
print(sbol3.get_namespace())
```

**Example output:**

```text
https://example.org/lab
None
```

## Other unwrapped rows

- sbol3: `document.data_path` is one row, also not wrapped; the other 117 rows are covered in section 2.
- screed: `hava_iter` and `read_hava_sequences` are two rows for the HAVA example format, `not-dna`.
- pyfaidx: `getmtime_fsspec` and `get_valid_filename` are two rows for file utilities, `not-dna`.
