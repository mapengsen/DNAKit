# sequence generation

Generate new DNA sequences through mutations, insertions, deletions, fragment rearrangements, random shuffling, or sequence recombination.


## 1) mutation generation {#1-ops-0051-mutation-generation}

- **Purpose:** Generate base substitutions with a per-base probability.
- **Method:** Independently select each base with mut_frac, then uniformly choose one of the other three bases.
- **API:** `dnakit.ops.evolution_generate(sequence[required], augmentations=("mutation",), seed/rng[choose one], mut_frac[optional])`.
- **Input:** A canonical DNASequence; probability in [0,1], default 0.05; provide either seed or rng.
- **Output:** Returns `EvolutionGenerationResult` with the generated `sequence` and mutation `steps`.
- **Example code:**

```python
from dnakit import DNASequence
from dnakit.ops import evolution_generate

seq = DNASequence("ACGTACGTACGT")
mutation_probability = 0.25
result = evolution_generate(
    seq,
    augmentations=("mutation",),
    seed=7,
    mut_frac=mutation_probability,  # Each base has an independent 25% mutation probability
)
print(result.sequence.symbols)
print(result.steps[0].length)  # Number of bases actually mutated
```

- **Example output:**

```text
ACATTCGACCTT
5
```

## 2) Insert generation {#2-ops-0052-insert-generation}

- **Purpose:** Insert random DNA segments with a per-position probability.
- **Method:** Select positions independently with insert_frac and insert random segments of length [insert_min,insert_max] after them.
- **API:** `dnakit.ops.evolution_generate(sequence[required], augmentations=("insertion",), seed/rng[choose one], insert_frac[optional], insert_min[optional], insert_max[optional])`.
- **Input:** A canonical DNASequence; probability in [0,1], positive insertion lengths and either seed or rng.
- **Output:** Returns `EvolutionGenerationResult` with the inserted `sequence` and insertion `steps`.
- **Example code:**

```python
from dnakit import DNASequence
from dnakit.ops import evolution_generate

seq = DNASequence("ACGTACGTACGT")
insertion_probability = 0.25

single_base = evolution_generate(
    seq,
    augmentations=("insertion",),
    seed=7,
    insert_frac=insertion_probability,  # Each position has an independent 25% insertion probability
    insert_min=1,
    insert_max=1,
)
segment = evolution_generate(
    seq,
    augmentations=("insertion",),
    seed=7,
    insert_frac=insertion_probability,
    insert_min=2,
    insert_max=4,
)
print(single_base.sequence.symbols, single_base.steps[0].length)
print(segment.sequence.symbols, segment.steps[0].length)
```

- **Example output:**

```text
ACGGTACTGTAACAGT 4
ACGAGACTTACAAACCGTTACAACCAGGT 17
```

## 3) delete generation {#3-ops-0053-delete-generation}

- **Purpose:** Delete bases with a per-base probability.
- **Method:** Independently select and remove each base with delete_frac.
- **API:** `dnakit.ops.evolution_generate(sequence[required], augmentations=("deletion",), seed/rng[choose one], delete_frac[optional])`.
- **Input:** A canonical DNASequence; probability in [0,1], default 0.05; provide either seed or rng.
- **Output:** Returns `EvolutionGenerationResult` with the deleted `sequence` and deletion `steps`.
- **Example code:**

```python
from dnakit import DNASequence
from dnakit.ops import evolution_generate

seq = DNASequence("ACGTACGTACGT")
deletion_probability = 0.25
result = evolution_generate(
    seq,
    augmentations=("deletion",),
    seed=7,
    delete_frac=deletion_probability,  # Each base has an independent 25% deletion probability
)
print(result.sequence.symbols)
print(result.steps[0].length)  # Number of bases actually deleted
```

- **Example output:**

```text
ACTCGCT
5
```

## 4) Evolutionary algorithm-like operation {#4-ops-0054-evolutionary-algorithm-like-operation}

- **Purpose:** Combine multiple random sequence transformations.
- **Method:** Sample operations without replacement and apply inversion, deletion, translocation, insertion, reverse complement and mutation in that order.
- **API:** `dnakit.ops.evolution_generate(sequence[required], augmentations[optional], max_augmentations[optional], hard_aug[optional], seed[optional], rng[optional], mut_frac[optional], insert_frac[optional], delete_frac[optional], insert_min/insert_max[optional], shift_min/shift_max[optional], invert_min/invert_max[optional], rc_prob[optional])`.
- **Input:** A canonical DNASequence, candidate operations and seed/rng; selected operations may leave the sequence unchanged.
- **Output:** Returns `EvolutionGenerationResult` with the generated `sequence` and operation `steps`.
- **Example code:**

```python
from dnakit import DNASequence
from dnakit.ops import evolution_generate

result = evolution_generate(
    DNASequence("AAGTCC"),
    augmentations=("reverse_complement",),
    rc_prob=1.0,
    seed=19,
)
print(result.sequence.symbols)
```

- **Example output:**

```text
GGACTT
```

## 5) Fragment rearrangement <span id="5"></span> {#5-ops-0055-fragment-rearrangement}

- **Purpose:** Randomly exchange, reverse-complement or duplicate segments.
- **Method:** Split at random breakpoints, then exchange segments, reverse-complement one or duplicate one according to operation.
- **API:** `dnakit.ops.rearrange_generate(sequence[required], operation[optional], segment_count[optional], seed[optional], rng[optional])`.
- **Input:** A canonical DNASequence, segment count, operation and seed/rng.
- **Output:** Returns a `RearrangementResult`; the example shows `sequence`, `breakpoints`, `permutation`.
- **Example code:**

```python
from dnakit import DNASequence
from dnakit.ops import rearrange_generate

seq = DNASequence("AAGTCC")
result = rearrange_generate(
    seq,
    operation="exchange",
    segment_count=3,
    seed=19,
)
print(result.sequence.symbols)
print(result.breakpoints, result.permutation)
```

- **Example output:**

```text
AGTCAC
(1, 5) (1, 0, 2)
```

## 6) k-mer randomly scrambled <span id="6-k-mer"></span> {#6-ops-0056-k-mer-randomly-scrambled}

- **Purpose:** Shuffle the sequence while preserving overlapping k-mer counts for the selected k.
- **Method:** Sample Eulerian paths in a de Bruijn multigraph to preserve overlapping k-mer counts exactly.
- **API:** `dnakit.ops.kmer_shuffle(sequence[required], k[optional], seed[optional], rng[optional], ensure_different[optional], max_attempts[optional])`.
- **Input:** A canonical DNASequence, k and seed/rng; by default a different sequence is required and impossible alternatives raise an error.
- **Output:** Returns a `KmerShuffleResult`; the example shows `sequence`.
- **Example code:**

```python
from collections import Counter

from dnakit import DNASequence
from dnakit.ops import kmer_shuffle

seq = DNASequence("AAGATCGATCGGATC")
result = kmer_shuffle(seq, k=3, seed=19)
generated = result.sequence.symbols
assert Counter(
    generated[index : index + 3] for index in range(len(generated) - 3 + 1)
) == Counter(seq.symbols[index : index + 3] for index in range(len(seq) - 3 + 1))
print(generated)
```

- **Example output:**

```text
AAGATCGGATCGATC
```

## 7) Multiple sequence recombination/crossover<span id="7-crossover"></span> {#7-ops-0057-multiple-sequence-recombinationcrossover}

- **Purpose:** Join two parent sequences at a single crossover boundary.
- **Method:** Take the first parent prefix and second parent suffix at the 0-based position boundary.
- **API:** `dnakit.ops.crossover(first[required], second[required], position[optional], seed[optional], rng[optional])`.
- **Input:** Two equal-length linear gap-free A/C/G/T sequences; provide seed/rng when position is omitted.
- **Output:** Returns a `CrossoverResult`; the example shows `sequence`.
- **Example code:**

```python
from dnakit import DNASequence
from dnakit.ops import crossover

first = DNASequence("AAAACCCC")
second = DNASequence("GGGGTTTT")
result = crossover(first, second, position=4)
print(result.sequence.symbols)  # AAAATTTT
```

- **Example output:**

```text
AAAATTTT
```

[FAQ: details](../../faq/tutorial_basics.md#generation-details)
