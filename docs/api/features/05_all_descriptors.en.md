# DNA descriptor

Calculate basic descriptors such as length, base composition, GC/AT, skewness, CpG, k-mer, entropy, complexity, repeats, and codons for DNA sequences. [FAQ: details](../../faq/tutorial_integrations.md#descriptors).

## 1) Fuzzy base statistics {#1-std-005-fuzzy-base-statistics}

- **Purpose:** Counts the number, proportion, and specific positions of `N` and other IUPAC ambiguous bases to assess sequence certainty and filter data by thresholds.
- **Method:** Identify ambiguous IUPAC symbols and report their counts, fractions and positions.
- **API:** `dnakit.normalize(raw[required], keep_ambiguous[optional], keep_u[optional], keep_other[optional], config[optional])`, `dnakit.validate(sequence[required], config[optional])`; the result field is `ambiguity` (`AmbiguityReport`).
- **Input:** Required IUPAC DNA; optional scaling denominator whether Gap is included, and alphabet and ambiguous base strategies.
- **Output:** Returns `NormalizationResult`; `ambiguity` contains ambiguous-base count, fraction, and counts by symbol.
- **Example code:**

```python
from dnakit import normalize

result = normalize("ANRY")
print(result.ambiguity.total_count)
print(result.ambiguity.fraction)
print([(item.symbol, item.count) for item in result.ambiguity.by_symbol])
```

- **Example output:**

```text
3
0.75
[('N', 1), ('R', 1), ('Y', 1)]
```

## 2) length characteristics {#2-desc-001-length-characteristics}

- **Purpose:** Calculate the base length of a DNA sequence and distinguish the effect of the number of symbols and Gap on the coordinate span for length screening, binning and window parameter settings.
- **Method:** Count sequence symbols separately from coordinate span including gaps.
- **API:** `dnakit.descriptors.length_features(value[required])`
- **Input:** Required `DNASequence` or `DNARecord`.
- **Output:** Returns a `LengthResult`; the example shows `canonical_base_count`.
- **Example code:**

```python
import dnakit
from dnakit.descriptors import length_features

seq = dnakit.normalize("ACGTACGT").sequence
assert seq is not None
result = length_features(seq)
print(len(seq))                     # 8
print(result.canonical_base_count)  # 8
```

- **Example output:**

```text
8
8
```

## 3) base composition {#3-desc-002-base-composition}

- **Purpose:** Count the number and proportion of A, C, G, and T in DNA sequences to form a basic composition vector for GC, skew, data distribution comparison and modeling.
- **Method:** Count A/C/G/T and use the selected ambiguity policy to determine denominators.
- **API:** `dnakit.descriptors.base_composition(value[required], ambiguity_policy[optional])`
- **Input:** Required `DNASequence` or `DNARecord`; optional `ambiguity_policy="error"|"ignore"`.
- **Output:** Returns a `CompositionResult`; the example shows `counts`, `fractions`.
- **Example code:**

```python
from dnakit import DNASequence
from dnakit.descriptors import base_composition

result = base_composition(DNASequence("AACGT"))
print(dict(result.counts), result.fractions["A"])
```

- **Example output:**

```text
{'A': 2, 'C': 1, 'G': 1, 'T': 1} 0.4
```

## 4) GC/AT features {#4-desc-003-gcat-features}

- **Purpose:** Calculate the number and ratio of GC and AT of DNA sequences, which is used to compare sequence composition, evaluate amplification or sequencing preferences, and serve as the basic input for other descriptors.
- **Method:** Sum G+C and A+T counts, then divide by the valid-base count.
- **API:** `dnakit.descriptors.gc_at_content(value[required], ambiguity_policy[optional])`
- **Input:** Required `DNASequence` or `DNARecord`; optional ambiguous-base handling policy.
- **Output:** Returns a `ContentResult`; the example shows `gc_fraction`, `at_fraction`.
- **Example code:**

```python
from dnakit import DNASequence
from dnakit.descriptors import gc_at_content

result = gc_at_content(DNASequence("AACG"))
print(result.gc_fraction, result.at_fraction)
```

- **Example output:**

```text
0.5 0.5
```

## 5) base skew {#5-desc-004-base-skew}

- **Purpose:** Calculate `(G-C)/(G+C)` and `(A-T)/(A+T)` respectively, quantify the degree of asymmetry of complementary bases on the current chain, and be used to observe local or overall composition bias.
- **Method:** Compute `(G-C)/(G+C)` and `(A-T)/(A+T)`.
- **API:** `dnakit.descriptors.base_skew(value[required], ambiguity_policy[optional])`
- **Input:** Required `DNASequence` or `DNARecord`; optional ambiguous-base handling policy.
- **Output:** Returns a `SkewResult`; the example shows `gc_skew`, `at_skew`.
- **Example code:**

```python
from dnakit import DNASequence
from dnakit.descriptors import base_skew

result = base_skew(DNASequence("AACG"))
print(result.gc_skew, result.at_skew)
```

- **Example output:**

```text
0.0 1.0
```

## 6) CpG characteristics {#6-desc-005-cpg-characteristics}

- **Purpose:** Count the number, density, and observed/expected values of CpG sites, and quantify the degree of CpG enrichment or deletion for CpG region screening and set comparison.
- **Method:** Count adjacent CpG pairs and derive density and observed/expected ratios from valid length and C/G counts.
- **API:** `dnakit.descriptors.cpg_features(value[required], ambiguity_policy[optional], cross_gaps[optional])`
- **Input:** Required `DNASequence` or `DNARecord`; optional `ambiguity_policy`, `cross_gaps`.
- **Output:** Returns a `CpGResult`; the example shows `cpg_count`, `density`, `observed_expected`.
- **Example code:**

```python
from dnakit import DNASequence
from dnakit.descriptors import cpg_features

result = cpg_features(DNASequence("ACGCGT"))
print(result.cpg_count, result.density, result.observed_expected)
```

- **Example output:**

```text
2 0.4 3.0
```

## 7) k-mer statistics {#7-desc-006-k-mer-statistics}

- **Purpose:** Enumerate k-mers of a specified length, return the number, frequency and existence, optionally merge reverse complementary k-mers for composition analysis, similarity calculation and feature modeling.
- **Method:** Scan k-base windows, optionally allowing overlap or merging reverse-complement words.
- **API:** `dnakit.descriptors.kmer_statistics(value[required], k[required], overlapping[optional], canonical[optional], ambiguity_policy[optional], cross_gaps[optional])`
- **Input:** Required `DNASequence`/`DNARecord` and `k`; optional overlap, canonical, ambiguous base and gap strategies.
- **Output:** Returns a `KmerResult`; the example shows `counts`, `denominator`.
- **Example code:**

```python
from dnakit import DNASequence
from dnakit.descriptors import kmer_statistics

result = kmer_statistics(DNASequence("ACGT"), 2)
print(dict(result.counts), result.denominator)
```

- **Example output:**

```text
{'AC': 1, 'CG': 1, 'GT': 1} 3
```

## 8) Sequence Entropy {#8-desc-007-sequence-entropy}

- **Purpose:** Calculate Shannon entropy based on the probability distribution of bases or k-mers, quantify the uniformity of distribution, and be used to identify sequences with low information content or single composition.
- **Method:** Compute `−Σ p log(p)` from base or k-mer frequencies.
- **API:** `dnakit.descriptors.shannon_entropy(value[required], unit[optional], k[optional], log_base[optional], ambiguity_policy[optional], cross_gaps[optional])`
- **Input:** Required `DNASequence`/`DNARecord`; optional `unit`, `k`, `log_base`, ambiguous-base and Gap policies.
- **Output:** Returns an `EntropyResult`; the example shows `entropy`.
- **Example code:**

```python
from dnakit import DNASequence
from dnakit.descriptors import shannon_entropy

result = shannon_entropy(DNASequence("ACGT"))
print(result.entropy)
```

- **Example output:**

```text
2.0
```

## 9) sequence complexity {#9-desc-008-sequence-complexity}

- **Purpose:** Compare the actual number of types of short sequences of different lengths with the theoretical number of possible types to obtain linguistic complexity, which is used to find sequences with more repetitions or insufficient pattern types.
- **Method:** Multiply observed/possible word-count ratios across word lengths.
- **API:** `dnakit.descriptors.linguistic_complexity(value[required], max_word_size[optional], ambiguity_policy[optional], cross_gaps[optional], max_observations[optional])`
- **Input:** Required `DNASequence`/`DNARecord`; optional maximum word length, ambiguous bases, Gap and workload cap.
- **Output:** Returns `ComplexityResult`; `score` is the sequence-complexity score.
- **Example code:**

```python
from dnakit import DNASequence
from dnakit.descriptors import linguistic_complexity

simple = linguistic_complexity(DNASequence("AAAAAAAA"), max_word_size=3)
diverse = linguistic_complexity(DNASequence("ACGTAGCT"), max_word_size=3)
print(simple.score < diverse.score)
```

- **Example output:**

```text
True
```

## 10) Homopolymer {#10-desc-009-homopolymer}

- **Purpose:** Search for consecutive repeats of the same base, return the position, base and length of each segment, and summarize the longest segment for screening sequencing, synthesis and amplification risks.
- **Method:** Scan runs of identical bases and retain those meeting the minimum length.
- **API:** `dnakit.descriptors.homopolymer_runs(value[required], min_run_length[optional], ambiguity_policy[optional], cross_gaps[optional])`
- **Input:** Required `DNASequence`/`DNARecord`; optional minimum continuous length, ambiguous base and Gap strategy.
- **Output:** Returns a `HomopolymerResult`; the example shows `longest_length`, `runs`.
- **Example code:**

```python
from dnakit import DNASequence
from dnakit.descriptors import homopolymer_runs

result = homopolymer_runs(DNASequence("AAACCGTTTT"), min_run_length=2)
print(result.longest_length, result.runs[-1].base)
```

- **Example output:**

```text
4 T
```

## 11) repetition ratio {#11-desc-010-repetition-ratio}

- **Purpose:** Find adjacent repeated sequence units and return repeated segments and coverage ratios, which are used to quantify highly repetitive regions and filter sequences that may affect alignment or synthesis.
- **Method:** Find exact tandem repeats and divide their union coverage by valid length.
- **API:** `dnakit.descriptors.exact_repeat_fraction(value[required], min_unit_length[optional], max_unit_length[optional], min_repeats[optional], ambiguity_policy[optional], cross_gaps[optional], max_comparisons[optional])`
- **Input:** Required `DNASequence`/`DNARecord`; optional repeat unit length, minimum number of times, Gap policy and comparison upper limit.
- **Output:** Returns an `ExactRepeatResult`; the example shows `repeat_fraction`, `runs`.
- **Example code:**

```python
from dnakit import DNASequence
from dnakit.descriptors import exact_repeat_fraction

result = exact_repeat_fraction(DNASequence("ATATATGC"), min_unit_length=2)
print(result.repeat_fraction, result.runs[0].unit)
```

- **Example output:**

```text
0.75 AT
```

## 12) window descriptor {#12-desc-011-window-descriptor}

- **Purpose:** Split DNA according to window length and step size, calculate GC, entropy and other descriptors window by window and retain positions, which is used to observe local changes in sequence features along coordinates.
- **Method:** Slide windows at the chosen step and calculate selected descriptors within each window.
- **API:** `dnakit.descriptors.window_descriptors(value[required], descriptors[required], window_size[required], step[optional], include_partial[optional], entropy_log_base[optional], ambiguity_policy[optional], cross_gaps[optional])`
- **Input:** Required sequence, descriptor list, and `window_size`; optional step size, `include_partial`, entropy base, ambiguous bases, and gap policy.
- **Output:** Returns a `WindowDescriptorResult`; the example shows `windows`.
- **Example code:**

```python
from dnakit import DNASequence
from dnakit.descriptors import window_descriptors

result = window_descriptors(
    DNASequence("ACGT"), ["gc", "entropy", "cpg"], window_size=2, step=2
)
print(result.windows[0].symbol_start, result.windows[0].values["gc_fraction"])
```

- **Example output:**

```text
0 0.5
```

## 13) codon statistics {#13-desc-012-codon-statistics}

- **Purpose:** Count codons and start and stop codons according to the specified reading frame, and return the count and position, which is used to analyze coding composition, reading frame and codon usage.
- **Method:** Count triplets in the selected frame and identify start/stop codons using the genetic code.
- **API:** `dnakit.descriptors.codon_statistics(value[required], frame[optional], genetic_code[optional], ambiguity_policy[optional], cross_gaps[optional])`
- **Input:** Required `DNASequence`/`DNARecord`; optional `frame=0|1|2`, `genetic_code=1`, ambiguous base and Gap strategies.
- **Output:** Returns a `CodonResult`; the example shows `codon_count`, `start_count`, `stop_count`.
- **Example code:**

```python
from dnakit import DNASequence
from dnakit.descriptors import codon_statistics

result = codon_statistics(DNASequence("ATGAAATAA"))
print(result.codon_count, result.start_count, result.stop_count)
```

- **Example output:**

```text
3 1 1
```

---


## 14) All descriptors (240 fields) {#all-descriptors}

- **Purpose:** Compute a fixed-order vector of composition, statistical and physicochemical descriptors.
- **Method:** Apply `descriptor_schema_v1`; unavailable fields are `None` with reasons.
- **API:** `dnakit.descriptors.all_descriptors(value, *, ambiguity_policy="ignore", conditions=None, dinucleotide_property_table=None)`
- **Input:** A `DNASequence` or `DNARecord`; the 60 dinucleotide-property fields require a user-supplied table. Calculation domains and default conditions are in [FAQ: details](../../faq/tutorial_integrations.md#descriptors).
- **Output:** Returns an `AllDescriptorsResult`; the example shows `values`, `unavailable_reasons`.
- **Example code:**

```python
from dnakit import DNASequence
from dnakit.descriptors import all_descriptors

result = all_descriptors(DNASequence("ACGT"))
print(len(result.values))
print(result.values["epsilon260_ss_m_inverse_cm_inverse"])
print(sum(name.startswith("diprodb_") and value is None for name, value in result.values.items()))
print(result.unavailable_reasons["diprodb_twist_mean"])
```

- **Example output:**

```text
240
40300.0
60
requires an explicit user-supplied DinucleotidePropertyTable; DNAKit bundles no DiProDB numerical values
```

<!-- BEGIN GENERATED API: sequence-descriptors -->

<span id="sequence-descriptors"></span>

## 15) Count a nucleotide or amino acid residue {#dnakit.databases.sequence_count_residues}

- **Purpose:** Count a nucleotide or amino acid residue.
- **Method:** Count the selected residue and report its fraction and 1-based positions.
- **API：** `dnakit.databases.sequence_count_residues(*, operation=None, sequence=None, uniprot_id=None, residue, verbose=True)`
- **Input:** Provide `residue`. Also provide sequence or uniprot_id; UniProt lookup requires network access.
- **Output:** Returns a mapping with `data` and related fields.
- **Example code:**

```python
import dnakit

result = dnakit.databases.sequence_count_residues(
    operation='count_residues',
    uniprot_id='P00533',
    residue='C',
    verbose=False,
)
print(result)
```

- **Example output:** The response has the following shape; records depend on the input and source.

Response-shape excerpt; `<…>` denotes type placeholders, not a captured query result.

```json
{
  "data": {
    "count": "<integer>",
    "sequence_length": "<integer>",
    "fraction": "<number>",
    "percent": "<number>",
    "positions_1based": [
      "<value>"
    ]
  }
}
```

[Full parameters and usage notes](../../faq/api/05_all_descriptors.md#dnakit.databases.sequence_count_residues-parameters)

## 16) Calculate GC content of a nucleotide sequence {#dnakit.databases.sequence_gc_content}

- **Purpose:** Calculate GC content of a nucleotide sequence.
- **Method:** Divide the G+C count by sequence length minus the N count.
- **API：** `dnakit.databases.sequence_gc_content(*, operation=None, sequence, verbose=True)`
- **Input:** Provide `sequence`.
- **Output:** Returns a mapping with `data` and related fields.
- **Example code:**

```python
import dnakit

result = dnakit.databases.sequence_gc_content(
    operation='gc_content',
    sequence='ATGCGATCGATCGATCGATCG',
    verbose=False,
)
print(result)
```

- **Example output:** The response has the following shape; records depend on the input and source.

Response-shape excerpt; `<…>` denotes type placeholders, not a captured query result.

```json
{
  "data": {
    "gc_fraction": "<number>",
    "gc_percent": "<number>",
    "length": "<integer>",
    "composition": "<object>"
  }
}
```

[Full parameters and usage notes](../../faq/api/05_all_descriptors.md#dnakit.databases.sequence_gc_content-parameters)

## 17) Summarize sequence length, composition and type-specific properties {#dnakit.databases.sequence_stats}

- **Purpose:** Summarize sequence length, composition and type-specific properties.
- **Method:** Count sequence length and composition, then calculate nucleotide- or protein-specific properties.
- **API：** `dnakit.databases.sequence_stats(*, operation=None, sequence=None, uniprot_id=None, verbose=True)`
- **Input:** Choose filters as shown in the example; omitted options use the operation defaults. Also provide sequence or uniprot_id; UniProt lookup requires network access.
- **Output:** Returns a mapping with `data` and related fields.
- **Example code:**

```python
import dnakit

result = dnakit.databases.sequence_stats(
    operation='stats',
    sequence='MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQVKVKALPDAQFEVVHSLAK',
    verbose=False,
)
print(result)
```

- **Example output:** The response has the following shape; records depend on the input and source.

Response-shape excerpt; `<…>` denotes type placeholders, not a captured query result.

```json
{
  "data": {
    "sequence_type": "<string>",
    "length": "<integer>",
    "composition": "<object>"
  }
}
```

[Full parameters and usage notes](../../faq/api/05_all_descriptors.md#dnakit.databases.sequence_stats-parameters)

<!-- END GENERATED API: sequence-descriptors -->
