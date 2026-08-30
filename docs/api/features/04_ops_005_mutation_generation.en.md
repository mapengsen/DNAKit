# OPS-005 sequence generation

Generate new DNA sequences through mutations, insertions, deletions, fragment rearrangements, random shuffling, or sequence recombination.

## 1. OPS-005.1 mutation generation

- **Function:** Independently determine whether a substitution mutation has occurred for each base in the sequence.
- **API:** `dnakit.ops.evolution_generate(sequence[required], augmentations=("mutation",), seed/rng[choose one], mut_frac[optional])`.
- **Probability:** `mut_frac` represents the probability of mutation of each base, the range is `0.0～1.0`, and the default value is `0.05`. The hit base will be replaced with one of the other three canonical bases with equal probability, so it will definitely change.
- **Supplementary:** `mutate()` is used for a single single base substitution at a specified position or a random position, and does not receive a probability parameter.
- **Sample code:**

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

- **Example results:**

```text
ACATTCGACCTT
5
```

## 2. OPS-005.2 Insert generation

- **Function:** Independently determine whether to insert random DNA after each base in the sequence.
- **API:** `dnakit.ops.evolution_generate(sequence[required], augmentations=("insertion",), seed/rng[choose one], insert_frac[optional], insert_min[optional], insert_max[optional])`.
- **Probability:** `insert_frac` represents the probability of insertion at each base position, the range is `0.0～1.0`, and the default value is `0.05`.
- **Insertion length:** `insert_min=1, insert_max=1`, only 1 random base will be inserted per hit; when other minimum and maximum values ​​are set, the fragment length will be randomly selected within the closed interval `[insert_min, insert_max]` for each hit.
- **Sample code:**

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

- **Example results:**

```text
ACGGTACTGTAACAGT 4
ACGAGACTTACAAACCGTTACAACCAGGT 17
```

If you just want to randomly select a position in the entire sequence and insert a continuous segment once, use `indel_generate(operation="insertion", min_length=..., max_length=...)`.

## 3. OPS-005.3 delete generation

- **Function:** Independently determine whether to delete the base for each base in the sequence.
- **API:** `dnakit.ops.evolution_generate(sequence[required], augmentations=("deletion",), seed/rng[choose one], delete_frac[optional])`.
- **Probability:** `delete_frac` represents the probability of deletion at each base position, the range is `0.0～1.0`, and the default value is `0.05`. After a hit, only one base at the current position is deleted.
- **Sample code:**

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

- **Example results:**

```text
ACTCGCT
5
```

If you only want to randomly select a contiguous segment in the entire sequence and delete it once, use `indel_generate(operation="deletion", min_length=..., max_length=...)`.

## 4. OPS-005.4 Evolutionary algorithm-like operation

`evolution_generate()` also supports partial reverse complementation, translocation, and entire reverse complementation. Multiple candidate operations are first sampled without replacement, and then executed in a fixed order: local reverse complementation, deletion, translocation, insertion, entire reverse complementation, and mutation. The base-by-base probabilities are for the current sequence when each operation is performed.

- **Full API:** `dnakit.ops.evolution_generate(sequence[required], augmentations[optional], max_augmentations[optional], hard_aug[optional], seed[optional], rng[optional], mut_frac[optional], insert_frac[optional], delete_frac[optional], insert_min/insert_max[optional], shift_min/shift_max[optional], invert_min/invert_max[optional], rc_prob[optional])`.
- **Result:** `EvolutionGenerationResult.sequence` is the generated sequence; `steps` saves whether each operation is effective and its actual impact length; the result also saves the random source, random state and `dnakit-evoaug-v3` algorithm version.
- **Sampling rules:** For `hard_aug=True`, exactly `max_augmentations` operations are selected for each sequence; for `False`, select randomly from 1 to the maximum value. Selecting a certain operation does not guarantee that the base will be hit. For example, if the probability is low and the sequence is short, no change may occur.

| Operation name | Function |
|----------------------|------------------------------------------------|
| `mutation` | Each base is independently replaced with `mut_frac` probability |
| `insertion` | After each base, random DNA is inserted independently with `insert_frac` probability |
| `deletion` | Each base is independently deleted with `delete_frac` probability |
| `translocation` | Random roll sequence, swap the order of fragments on both sides of the breakpoint |
| `inversion` | Replace contiguous segments with their reverse complements |
| `reverse_complement` | Press `rc_prob` to replace the entire sequence with the reverse complement |

## 5. OPS-005.5 Fragment rearrangement <span id="5"></span>

`rearrange_generate()` First use random breakpoints to cut the sequence into `segment_count` non-empty segments, and then perform a rearrangement:

| `operation` | Function |
| --------------- | ----------------------------------------------- |
| `exchange` | Randomly swap the order of fragments and keep the length unchanged |
| `inversion` | Select a segment and replace it with its reverse complement, keeping the length unchanged |
| `duplication` | Select a segment and make a copy after it, increasing the length |

- **API:** `dnakit.ops.rearrange_generate(sequence[required], operation[optional], segment_count[optional], seed[optional], rng[optional])`.
- **Result:** `RearrangementResult` Save `breakpoints`, swap `permutation` or `selected_segment` for easy review of the build process.

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

## 6. OPS-005.6 k-mer randomly scrambled <span id="6-k-mer"></span>

`kmer_shuffle()` Instead of directly scrambling bases independently, Euler paths are randomly sampled on the de Bruijn multigraph of the original sequence. Thus the count of overlapping k-mers is maintained exactly; for example, `k=2` maintains the count of doublets, and `k=3` maintains the count of triplets. By default `ensure_different=True`, an attempt will be made to generate a sequence that is different from the original sequence; if the counting structure has only a unique reconstruction result, an explicit error will be thrown.

- **API:** `dnakit.ops.kmer_shuffle(sequence[required], k[optional], seed[optional], rng[optional], ensure_different[optional], max_attempts[optional])`.
- **Results:** `KmerShuffleResult` Saves `k`, raw k-mer count, number of attempts and random state.

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

This implementation refers to the k-mer-preserving shuffled Euler path method (see [uShuffle paper ](https://pmc.ncbi.nlm.nih.gov/articles/PMC2375906/)); it only guarantees the overlap count of the specified `k`, and does not mean to preserve higher-order k-mers, functionality, expression volume, or experimental properties.

## 7. OPS-005.7 Multiple sequence recombination/crossover<span id="7-crossover"></span>

`crossover()` Receives two DNA sequences of equal length and a single point boundary, and generates the children of "the first sequence prefix + the second sequence suffix". The boundary uses 0-based half-open splitting positions; when `position` is omitted, `seed` or `random.Random` must be passed, and a boundary is randomly selected that preserves the bases on both sides.

- **API:** `dnakit.ops.crossover(first[required], second[required], position[optional], seed[optional], rng[optional])`.
- **Result:** `CrossoverResult.sequence` is the offspring; `position`, the two parent lengths and the random state are saved in the result.
- **Example:**

```python
from dnakit import DNASequence
from dnakit.ops import crossover

first = DNASequence("AAAACCCC")
second = DNASequence("GGGGTTTT")
result = crossover(first, second, position=4)
print(result.sequence.symbols)  # AAAATTTT
```

The API is currently a one-point intersection of two sequences of equal length; inputs of different lengths, circular shapes, containing gaps, and containing IUPAC ambiguous bases will be rejected.
