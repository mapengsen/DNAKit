# Legality check

Checks a DNA sequence, record, or data set for compliance with rules such as alphabet, length, gap, quality, and metadata, and returns a structured list of questions.

Whenever an error is found, the `is_valid` of the returned result is `False`; all issues are still stored in the structured `issues`.

## 1) Validate sequences and datasets

- **Purpose:** Check alphabet, illegal symbols, empty sequences, length, gap, ambiguous base ratio, record metadata, PHRED quality, empty set and duplicate ID in one go, and return an itemized issue report to filter out or correct invalid data before analysis.
- **Method:** Evaluate each configured rule and collect issues; errors set is_valid=False while warnings alone do not.

- **API:** `dnakit.validate(value[required], config[optional])`; single input returns `ValidationReport`, collective input returns `DatasetValidationReport`. Configuration type is `dnakit.ValidationConfig` or `dnakit.DatasetValidationConfig`.
- **Input:** Ordinary users can directly input `DNA`. When there is only one record, use `ValidationConfig` and return a single report; when it contains multiple records, you can use `DatasetValidationConfig`, or you can directly pass `ValidationConfig` as the rule for each record. Old core objects are still compatible.
- **Output:** A single input returns `ValidationReport`; a collection returns `DatasetValidationReport`. Both expose validity and issues.
- **Example code:**

```python
from dnakit import DNA, ValidationConfig, validate

single = validate(
    DNA("ACGN", alphabet="iupac"),
    config=ValidationConfig(alphabet="strict"),
)
print(single.is_valid)
print([issue.code for issue in single.issues])

records = DNA(
    [
        {"sequence": "ACGT", "id": "valid"},
        {"sequence": "AC", "id": "short"},
    ]
)
collection = validate(records, config=ValidationConfig(min_length=3))
print(collection.is_valid)
print([issue.code for issue in collection.issues])
```

- **Example output:**

```text
False
['STD_INVALID_SYMBOL']
False
['STD_INVALID_RECORDS']
```

[FAQ: details](../../faq/tutorial_basics.md#validation-cli)
