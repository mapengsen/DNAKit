# normalize

Unify the handling of case, whitespace, `U/T` and illegal characters in DNA sequences, resulting in standardized sequences.

## 1) STD-001 character standardization

- **Function:** Organize the original text into a unified DNA sequence format, handle case, blank, `U`, ambiguous characters and illegal characters according to the settings, and record the specific modifications for consistent input in subsequent calculations.
- **API:** `dnakit.normalize(raw[required], keep_ambiguous[optional], keep_u[optional], keep_other[optional], config[optional])`; `config` uses `dnakit.NormalizationConfig`.
- **Input:** Required raw string, UTF-8 `bytes`, `DNASequence`, or iterable consisting of string and `Gap`.
- **Core args:**
  - `keep_ambiguous=True`: retain `RYSWKMBDHVN` ambiguous base; delete when set to `False`.
  - `keep_u=False`: Delete `U` by default; retain when set to `True`.
  - `keep_other=False`: Removes characters except `A/C/G/T`, IUPAC ambiguous bases, and `U` by default; retained when set to `True`.
- **Default behavior:** Keep ambiguous bases, remove `U` and other non-DNA characters.
- **Sample code:**

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

- **Example results:**

```text
ACN
iupac
['delete_other', 'delete_u', 'delete_other']
```
