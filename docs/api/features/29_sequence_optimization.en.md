# Sequence optimization

## 1) Optimize DNA under sequence constraints {#bio-45}

- **Purpose:** Modify DNA to satisfy rules, optionally preserving its protein translation.
- **Method:** Satisfy hard constraints before weighted objectives; this example removes `GCT` through synonymous codon changes.
- **API:** `dnakit.bioinformatics.optimization.specification(kind, region=None, weight=1, **parameters)`; `dnakit.bioinformatics.optimization.optimize_sequence(sequence, constraints=(), objectives=(), preserve_translation=False, seed=0)`
- **Input:** An A/C/G/T sequence and rules; length is fixed; see [GC, pattern, structure and codon rules](../../faq/bioinformatics.md#optimization).
- **Output:** Returns `OptimizationResult` with the optimized sequence, constraint status, and edits.
- **Example code:**

```python
from dnakit.bioinformatics.optimization import optimize_sequence, specification

result = optimize_sequence(
    "GCTGCT", preserve_translation=True,
    constraints=[specification("avoid_pattern", pattern="GCT")], seed=0,
)
print(result.sequence, result.constraints_satisfied)
print(result.edits)
```

- **Example output:**

```text
GCAGCA True
((2, 'T', 'A'), (5, 'T', 'A'))
```

[Dependencies and implementation notes](../../faq/bioinformatics.md#optimization).

## 2) Define custom constraints and objectives {#bio-46}

- **Purpose:** Add your own evaluation logic to sequence optimization.
- **Method:** Return `SpecEvaluation`; constraints minimize nonnegative `violation`, while objectives maximize `utility`.
- **API:** `dnakit.bioinformatics.optimization.Specification(name, evaluator, weight=1)`; `dnakit.bioinformatics.optimization.SpecEvaluation(violation=0, utility=0, locations=(), message='')`
- **Input:** A Python callable taking `(sequence, circular)` and returning finite evaluation values.
- **Output:** Returns `OptimizationResult` with the optimized sequence and constraint status.
- **Example code:**

```python
from dnakit.bioinformatics.optimization import Specification, SpecEvaluation, optimize_sequence

rule = Specification("starts_with_A",
                     lambda sequence, circular: SpecEvaluation(violation=float(sequence[0] != "A")))
result = optimize_sequence("CCGT", constraints=[rule])
print(result.sequence, result.constraints_satisfied)
```

- **Example output:**

```text
ACGT True
```

[Dependencies and implementation notes](../../faq/bioinformatics.md#optimization).

## 3) Read optimization rules from GenBank annotations {#bio-47}

- **Purpose:** Convert rule annotations in construct files into constraints and objectives.
- **Method:** Parse JSON in feature qualifiers `dnakit_constraint` / `dnakit_objective` and use the feature region.
- **API:** `dnakit.bioinformatics.optimization.specifications_from_record(record)`
- **Input:** A Biopython SeqRecord with exact, forward, simple feature locations; requires `bio-formats`.
- **Output:** Returns tuples of constraint and objective rules for sequence optimization.
- **Example code:**

```python
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, SimpleLocation
from dnakit.bioinformatics.optimization import specifications_from_record

record = SeqRecord(Seq("ACGTACGT"), id="construct")
record.features = [SeqFeature(SimpleLocation(0, 8, strand=1),
    qualifiers={"dnakit_constraint": ['{"kind":"gc","minimum":0.4,"maximum":0.6}']})]
constraints, objectives = specifications_from_record(record)
print(len(constraints), len(objectives), constraints[0].evaluate(str(record.seq)).violation)
```

- **Example output:**

```text
1 0 0.0
```

[Dependencies and implementation notes](../../faq/bioinformatics.md#optimization).

## 4) Constrain sequence design with BLAST hits {#bio-48}

- **Purpose:** Detect homologous HSPs against a chosen database during optimization.
- **Method:** Run BLAST and score violations by identity and hit-length thresholds, caching by sequence.
- **API:** `dnakit.bioinformatics.optimization.blast_match_specification(database, minimum_identity=0.8, minimum_length=20, evalue=1e-5)`
- **Input:** An existing local nucleotide BLAST database prefix and blastn; no hits do not prove the absence of all homology.
- **Output:** Returns a BLAST-match `Specification`; evaluation gives violation severity and hit intervals.
- **Example code:**

```python
from dnakit.bioinformatics.optimization import blast_match_specification

rule = blast_match_specification("local_db/database", minimum_identity=0.9, minimum_length=20)
evaluation = rule.evaluate("ACGTACGTACGTACGTACGTACGT")
print(evaluation.violation, evaluation.locations)
```

- **Example output:**

Returns `SpecEvaluation` with a violation score and hit locations; actual results depend on the database. This example was not executed.

[Dependencies and implementation notes](../../faq/bioinformatics.md#optimization).
