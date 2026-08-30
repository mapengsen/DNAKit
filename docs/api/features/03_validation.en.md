# Legality check

Checks a DNA sequence, record, or data set for compliance with rules such as alphabet, length, gap, quality, and metadata, and returns a structured list of questions.

Whenever an error is found, the `is_valid` of the returned result is `False`; all issues are still stored in the structured `issues`.

- **Use:** Check alphabet, illegal symbols, empty sequences, length, gap, ambiguous base ratio, record metadata, PHRED quality, empty set and duplicate ID in one go, and return an itemized issue report to filter out or correct invalid data before analysis.

- **API:** `dnakit.validate(value[required], config[optional])`; single input returns `ValidationReport`, collective input returns `DatasetValidationReport`. Configuration type is `dnakit.ValidationConfig` or `dnakit.DatasetValidationConfig`.
- **Input:** Ordinary users can directly input `DNA`. When there is only one record, use `ValidationConfig` and return a single report; when it contains multiple records, you can use `DatasetValidationConfig`, or you can directly pass `ValidationConfig` as the rule for each record. Old core objects are still compatible.
- **Check results:** `is_valid=False` for a single result indicates that at least one rule produced an error; `is_valid=False` for a collection result indicates that the collection is empty, any record is illegal, or a duplicate ID is found when unique ID checking is enabled. Warnings alone do not render the result illegal.
- **Sample code:**

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

- **Example results:**

```text
False
['STD_INVALID_SYMBOL']
False
['STD_INVALID_RECORDS']
```

- **Command line:** The original text command will still be standardized first, and then the unified verification logic will be called:

```bash
dnakit validate ACGT --sequence-length 4
```
