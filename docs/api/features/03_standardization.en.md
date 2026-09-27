# normalize

Unify the handling of case, whitespace, `U/T` and illegal characters in DNA sequences, resulting in standardized sequences.

## 1) character standardization {#1-std-001-character-standardization}

- **Purpose:** Organize the original text into a unified DNA sequence format, handle case, blank, `U`, ambiguous characters and illegal characters according to the settings, and record the specific modifications for consistent input in subsequent calculations.
- **Method:** Remove whitespace, normalize case and retain or remove U, ambiguous bases and other symbols according to the configuration; record each change.
- **API:** `dnakit.normalize(raw[required], keep_ambiguous[optional], keep_u[optional], keep_other[optional], config[optional])`; `config` uses `dnakit.NormalizationConfig`.
- **Input:** Required raw string, UTF-8 `bytes`, `DNASequence`, or iterable consisting of string and `Gap`. Defaults: `keep_ambiguous=True`, `keep_u=False`, `keep_other=False`.
- **Output:** Returns a `NormalizationResult`; the example shows `sequence`, `changes`.
- **Example code:**

```python
from dnakit import normalize

result = normalize(
    " acn-uX\n",
    keep_ambiguous=True,
    keep_u=False,
    keep_other=False,
)
print(result.sequence.symbols)
print(result.sequence.alphabet.value)
print(
    [
        change.operation
        for change in result.changes
        if change.operation.startswith("delete_")
    ]
)
```

- **Example output:**

```text
ACN
iupac
['delete_other', 'delete_u', 'delete_other']
```

[FAQ: details](../../faq/tutorial_basics.md#normalization-options)
