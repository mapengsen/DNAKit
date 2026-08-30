# DNA descriptor

Calculate basic descriptors such as length, base composition, GC/AT, skewness, CpG, k-mer, entropy, complexity, repeats, and codons for DNA sequences.

The result object can be converted into a serializable dictionary via `to_dict()`. Where the interface provides `cross_gaps`, `True` is only allowed to span spanable Gap; the position of `Gap(crossable=False)` explicitly set is always a hard boundary. Descriptors that accept circular sequences still use the current origin as the scan boundary, and do not automatically complement adjacent pairs or windows from the end to the beginning; when cross-origin semantics are required, the sequence should be rotated or explicitly linearized first.

A total of 13 basic functions that can be called individually are listed below: `STD-005` and `DESC-001`–`DESC-012`; when you need to output a more complete fixed feature vector at one time, please continue to view [ "All descriptor calculations" ](#all-descriptors) later in this page.

## 1) STD-005 Fuzzy base statistics

- **Role:** Counts the number, proportion, and specific positions of `N` and other IUPAC ambiguous bases to assess sequence certainty and filter data by thresholds.
- **API:** `dnakit.normalize(raw[required], keep_ambiguous[optional], keep_u[optional], keep_other[optional], config[optional])`, `dnakit.validate(sequence[required], config[optional])`; the result field is `ambiguity` (`AmbiguityReport`).
- **Input:** Required IUPAC DNA; optional scaling denominator whether Gap is included, and alphabet and ambiguous base strategies.
- **Sample code:**

```python
from dnakit import normalize

result = normalize("ANRY")
print(result.ambiguity.total_count)
print(result.ambiguity.fraction)
print([(item.symbol, item.count) for item in result.ambiguity.by_symbol])
```

- **Example results:**

```text
3
0.75
[('N', 1), ('R', 1), ('Y', 1)]
```

## 2) DESC-001 length characteristics

- **Function:** Calculate the base length of a DNA sequence and distinguish the effect of the number of symbols and Gap on the coordinate span for length screening, binning and window parameter settings.
- **API:** `dnakit.descriptors.length_features(value[required])`
- **Input:** Required `DNASequence` or `DNARecord`.
- **Sample code:**

```python
import dnakit
from dnakit.descriptors import length_features

seq = dnakit.normalize("ACGTACGT").sequence
assert seq is not None
result = length_features(seq)
print(len(seq))                     # 8
print(result.canonical_base_count)  # 8
```

- **Example results:**

```text
8
8
```

## 3) DESC-002 base composition

- **Function:** Count the number and proportion of A, C, G, and T in DNA sequences to form a basic composition vector for GC, skew, data distribution comparison and modeling.
- **API:** `dnakit.descriptors.base_composition(value[required], ambiguity_policy[optional])`
- **Input:** Required `DNASequence` or `DNARecord`; optional `ambiguity_policy="error"|"ignore"`.
- **Sample code:**

```python
from dnakit import DNASequence
from dnakit.descriptors import base_composition

result = base_composition(DNASequence("AACGT"))
print(dict(result.counts), result.fractions["A"])
```

- **Example results:**

```text
{'A': 2, 'C': 1, 'G': 1, 'T': 1} 0.4
```

## 4) DESC-003 GC/AT features

- **Function:** Calculate the number and ratio of GC and AT of DNA sequences, which is used to compare sequence composition, evaluate amplification or sequencing preferences, and serve as the basic input for other descriptors.
- **API:** `dnakit.descriptors.gc_at_content(value[required], ambiguity_policy[optional])`
- **Input:** Required `DNASequence` or `DNARecord`; optional fuzzy base handling strategy.
- **Sample code:**

```python
from dnakit import DNASequence
from dnakit.descriptors import gc_at_content

result = gc_at_content(DNASequence("AACG"))
print(result.gc_fraction, result.at_fraction)
```

- **Example results:**

```text
0.5 0.5
```

## 5) DESC-004 base skew

- **Function:** Calculate `(G-C)/(G+C)` and `(A-T)/(A+T)` respectively, quantify the degree of asymmetry of complementary bases on the current chain, and be used to observe local or overall composition bias.
- **API:** `dnakit.descriptors.base_skew(value[required], ambiguity_policy[optional])`
- **Input:** Required `DNASequence` or `DNARecord`; optional fuzzy base handling strategy.
- **Sample code:**

```python
from dnakit import DNASequence
from dnakit.descriptors import base_skew

result = base_skew(DNASequence("AACG"))
print(result.gc_skew, result.at_skew)
```

- **Example results:**

```text
0.0 1.0
```

## 6) DESC-005 CpG characteristics

- **Function:** Count the number, density, and observed/expected values of CpG sites, and quantify the degree of CpG enrichment or deletion for CpG region screening and set comparison.
- **API:** `dnakit.descriptors.cpg_features(value[required], ambiguity_policy[optional], cross_gaps[optional])`
- **Input:** Required `DNASequence` or `DNARecord`; optional `ambiguity_policy`, `cross_gaps`.
- **Sample code:**

```python
from dnakit import DNASequence
from dnakit.descriptors import cpg_features

result = cpg_features(DNASequence("ACGCGT"))
print(result.cpg_count, result.density, result.observed_expected)
```

- **Example results:**

```text
2 0.4 3.0
```

## 7) DESC-006 k-mer statistics

- **Function:** Enumerate k-mers of a specified length, return the number, frequency and existence, optionally merge reverse complementary k-mers for composition analysis, similarity calculation and feature modeling.
- **API:** `dnakit.descriptors.kmer_statistics(value[required], k[required], overlapping[optional], canonical[optional], ambiguity_policy[optional], cross_gaps[optional])`
- **Input:** Required `DNASequence`/`DNARecord` and `k`; optional overlap, canonical, ambiguous base and gap strategies.
- **Sample code:**

```python
from dnakit import DNASequence
from dnakit.descriptors import kmer_statistics

result = kmer_statistics(DNASequence("ACGT"), 2)
print(dict(result.counts), result.denominator)
```

- **Example results:**

```text
{'AC': 1, 'CG': 1, 'GT': 1} 3
```

## 8) DESC-007 Sequence Entropy

- **Function:** Calculate Shannon entropy based on the probability distribution of bases or k-mers, quantify the uniformity of distribution, and be used to identify sequences with low information content or single composition.
- **API:** `dnakit.descriptors.shannon_entropy(value[required], unit[optional], k[optional], log_base[optional], ambiguity_policy[optional], cross_gaps[optional])`
- **Input:** Required `DNASequence`/`DNARecord`; optional `unit`, `k`, `log_base`, fuzzy base and Gap strategy.
- **Sample code:**

```python
from dnakit import DNASequence
from dnakit.descriptors import shannon_entropy

result = shannon_entropy(DNASequence("ACGT"))
print(result.entropy)
```

- **Example results:**

```text
2.0
```

## 9) DESC-008 sequence complexity

- **Function:** Compare the actual number of types of short sequences of different lengths with the theoretical number of possible types to obtain linguistic complexity, which is used to find sequences with more repetitions or insufficient pattern types.
- **API:** `dnakit.descriptors.linguistic_complexity(value[required], max_word_size[optional], ambiguity_policy[optional], cross_gaps[optional], max_observations[optional])`
- **Input:** Required `DNASequence`/`DNARecord`; optional maximum word length, fuzzy bases, Gap and workload cap.
- **Sample code:**

```python
from dnakit import DNASequence
from dnakit.descriptors import linguistic_complexity

simple = linguistic_complexity(DNASequence("AAAAAAAA"), max_word_size=3)
diverse = linguistic_complexity(DNASequence("ACGTAGCT"), max_word_size=3)
print(simple.score < diverse.score)
```

- **Example results:**

```text
True
```

## 10) DESC-009 Homopolymer

- **Function:** Search for consecutive repeats of the same base, return the position, base and length of each segment, and summarize the longest segment for screening sequencing, synthesis and amplification risks.
- **API:** `dnakit.descriptors.homopolymer_runs(value[required], min_run_length[optional], ambiguity_policy[optional], cross_gaps[optional])`
- **Input:** Required `DNASequence`/`DNARecord`; optional minimum continuous length, ambiguous base and Gap strategy.
- **Sample code:**

```python
from dnakit import DNASequence
from dnakit.descriptors import homopolymer_runs

result = homopolymer_runs(DNASequence("AAACCGTTTT"), min_run_length=2)
print(result.longest_length, result.runs[-1].base)
```

- **Example results:**

```text
4 T
```

## 11) DESC-010 repetition ratio

- **Function:** Find adjacent repeated sequence units and return repeated segments and coverage ratios, which are used to quantify highly repetitive regions and filter sequences that may affect alignment or synthesis.
- **API:** `dnakit.descriptors.exact_repeat_fraction(value[required], min_unit_length[optional], max_unit_length[optional], min_repeats[optional], ambiguity_policy[optional], cross_gaps[optional], max_comparisons[optional])`
- **Input:** Required `DNASequence`/`DNARecord`; optional repeat unit length, minimum number of times, Gap policy and comparison upper limit.
- **Sample code:**

```python
from dnakit import DNASequence
from dnakit.descriptors import exact_repeat_fraction

result = exact_repeat_fraction(DNASequence("ATATATGC"), min_unit_length=2)
print(result.repeat_fraction, result.runs[0].unit)
```

- **Example results:**

```text
0.75 AT
```

## 12) DESC-011 window descriptor

- **Function:** Split DNA according to window length and step size, calculate GC, entropy and other descriptors window by window and retain positions, which is used to observe local changes in sequence features along coordinates.
- **API:** `dnakit.descriptors.window_descriptors(value[required], descriptors[required], window_size[required], step[optional], include_partial[optional], entropy_log_base[optional], ambiguity_policy[optional], cross_gaps[optional])`
- **Input:** Required sequence, descriptor list, and `window_size`; optional step size, `include_partial`, entropy base, fuzzy bases, and gap policy.
- **Sample code:**

```python
from dnakit import DNASequence
from dnakit.descriptors import window_descriptors

result = window_descriptors(
    DNASequence("ACGT"), ["gc", "entropy", "cpg"], window_size=2, step=2
)
print(result.windows[0].symbol_start, result.windows[0].values["gc_fraction"])
```

- **Example results:**

```text
0 0.5
```

## 13) DESC-012 codon statistics

- **Function:** Count codons and start and stop codons according to the specified reading frame, and return the count and position, which is used to analyze coding composition, reading frame and codon usage.
- **API:** `dnakit.descriptors.codon_statistics(value[required], frame[optional], genetic_code[optional], ambiguity_policy[optional], cross_gaps[optional])`
- **Input:** Required `DNASequence`/`DNARecord`; optional `frame=0|1|2`, `genetic_code=1`, ambiguous base and Gap strategies.
- **Sample code:**

```python
from dnakit import DNASequence
from dnakit.descriptors import codon_statistics

result = codon_statistics(DNASequence("ATGAAATAA"))
print(result.codon_count, result.start_count, result.stop_count)
```

- **Example results:**

```text
3 1 1
```

---



## 14) All descriptor calculations (240 items)

**Full 240 descriptors**

Compute the compositional, statistical, physicochemical, and dinucleotide descriptors of DNA sequences at once to obtain 240 feature results in a fixed order.

`all_descriptors()` Use the fixed version `descriptor_schema_v1` and return the non-computable reason, calculation conditions and source information at the same time. The order of the fields will not change with the input; the first 180 items are calculated according to the input and their respective applicable fields. DNAKit does not have a built-in dinucleotide value table, so the last 60 items default to `None`.

<span id="1"></span>**Shortest usage**

```python
from dnakit import DNASequence
from dnakit.descriptors import all_descriptors

result = all_descriptors(DNASequence("ACGT"))
print(len(result.values))
print(result.values["epsilon260_ss_m_inverse_cm_inverse"])
print(sum(name.startswith("diprodb_") and value is None for name, value in result.values.items()))
print(result.unavailable_reasons["diprodb_twist_mean"])
```

```text
240
40300.0
60
requires an explicit user-supplied DinucleotidePropertyTable; DNAKit bundles no DiProDB numerical values
```

The CLI outputs the same complete set of fields by default; if you need the old four sets of condensed results, you can use `--compact` explicitly:

```bash
dnakit describe ACGT
dnakit describe ACGT --compact
```

<span id="2-api"></span>**API**

- `dnakit.descriptors.all_descriptors(value[required], ambiguity_policy[optional], conditions[optional], dinucleotide_property_table[optional])`: Returns all fields; the last 60 are calculated only if the table is explicitly provided.
- `dnakit.descriptors.load_dinucleotide_property_table(path[required])`: Bounded read and strict validation of user JSON table, record file SHA-256.
- `dnakit.descriptors.descriptor_schema_v1()`: Returns 240 immutable `DescriptorField`, each containing `index`, `name`, `category`, `unit`, `formula`, and `source`.
- `dnakit.descriptors.DESCRIPTOR_SCHEMA_V1`: constant form of the same fixed schema.
- `dnakit.descriptors.DESCRIPTOR_NAMES_V1`: Contains only ordered field names.

Result object structure:

| Properties | Meaning |
| ----------------------- | --------------------------------------------------------------- |
| `schema_version` | Fixed to `descriptor_schema_v1` |
| `sequence_id` | Keep record ID if input is `DNARecord`, otherwise `None` |
| `values` | 240 values strictly arranged by schema |
| `unavailable_reasons` | Contains only fields with value `None`, and each field has exactly one reason |
| `conditions` | IUPAC, Gap, k-mer, ORF, replicate, thermodynamic and parameter table conditions |
| `provenance` | DNAKit version, implementation type, and user table declaration name, version, origin, and SHA-256 |

<span id="3"></span>**Unified calculation caliber**

- Only A/C/G/T counts in the canonical denominator; defaults to `ambiguity_policy="ignore"`, can also be set to `"error"`.
- k-mers, CpG/GpC, user dinucleotide tables, and other neighbor word statistics all allow overlap, but do not span IUPAC ambiguous symbols or explicit gaps.
- All scales returning `None` when the denominator is 0 are not replaced with fake `0`; the reason is written to `unavailable_reasons`.
- User table statistics `sd` use the population standard deviation; the table must provide a specified set of 15 attributes and a finite number of 16 DNA dinucleotides per set.
- frame 0 uses NCBI standard genetic code 1; the complete ORF uses `ATG` to start and `TAA/TAG/TGA` to stop, and scans the forward and reverse strands for a total of six reading frames.
- Thermodynamic default conditions are 37 °C, Na⁺ 0.05 M, chain concentration 250 nM; Wallace only 2–13 nt, SantaLucia NN only 2–60 nt.
- Molecular weight and ε260 are only for linear, no-gap, no-modification, A/C/G/T DNA; ε260 is a theoretical coefficient, not experimental A260.
- Currently LZ76 calculates sequences that are linear, have no gaps, no ambiguous symbols and do not exceed 10,000 nt; longer sequences return `None` to avoid runaway quadratic time complexity.

<span id="4"></span>**Category and field range**

| Serial number | Category | Number of fields |
| -------- | ---------------------------------------- | -----: |
| 1–12 | Base Length and Data Quality | 12 |
| 13–28 | Chemical grouping of bases | 16 |
| 29–112 | 1/2/3-mer frequency | 84 |
| 113–128 | skew, CpG/GpC and Chargaff metrics | 16 |
| 129–148 | Entropy, complexity, homopolymers and repetition | 20 |
| 149–164 | Codons and six reading frame ORF | 16 |
| 165–180 | Molecular weight, ε260, Tm and NN Thermodynamics | 16 |
| 181–240 | 15 sets of user dinucleotide parameters × mean/sd/min/max | 60 |

<span id="5-240"></span>**240 complete field table**

| # | Field | Category | Unit | Formula | Source |
| --: | -------------------------------------------------- | ---------------------------------- | --------------- | ----------------------------------------------------------------------------------------------- | ------------------------------------------------------------------------------- |
| 1 | `symbol_length` | `basic` | `nt` | `number of nucleotide symbols; explicit gaps excluded` | DNAKit descriptor_schema_v1 |
| 2 | `coordinate_span` | `basic` | `nt` | `symbol_length + sum(known gap lengths); undefined with an unknown gap` | DNAKit descriptor_schema_v1 |
| 3 | `canonical_base_count` | `basic` | `count` | `count(A,C,G,T)` | DNAKit descriptor_schema_v1 |
| 4 | `ambiguity_symbol_count` | `basic` | `count` | `symbol_length - canonical_base_count` | DNAKit descriptor_schema_v1 |
|   5 | `gap_object_count`                         | `basic`                 | `count`       | `number of explicit Gap objects`                                                           | DNAKit descriptor_schema_v1                                        |
|   6 | `known_gap_nt`                             | `basic`                 | `nt`          | `sum(length of known explicit gaps)`                                                       | DNAKit descriptor_schema_v1                                        |
|   7 | `unknown_gap_count`                        | `basic`                 | `count`       | `number of explicit gaps with unknown length`                                              | DNAKit descriptor_schema_v1                                        |
|   8 | `canonical_symbol_fraction`                | `basic`                 | `fraction`    | `canonical_base_count / symbol_length`                                                     | DNAKit descriptor_schema_v1                                        |
|   9 | `ambiguity_symbol_fraction`                | `basic`                 | `fraction`    | `ambiguity_symbol_count / symbol_length`                                                   | DNAKit descriptor_schema_v1                                        |
|  10 | `known_gap_fraction`                       | `basic`                 | `fraction`    | `known_gap_nt / coordinate_span`                                                           | DNAKit descriptor_schema_v1                                        |
|  11 | `canonical_run_count`                      | `basic`                 | `count`       | `number of uninterrupted A/C/G/T runs split by ambiguity or gaps`                          | DNAKit descriptor_schema_v1                                        |
|  12 | `longest_canonical_run_nt`                 | `basic`                 | `nt`          | `maximum uninterrupted A/C/G/T run length`                                                 | DNAKit descriptor_schema_v1                                        |
|  13 | `purine_count`                             | `composition`           | `count`       | `count(A)+count(G)`                                                                        | DNAKit descriptor_schema_v1                                        |
|  14 | `purine_fraction`                          | `composition`           | `fraction`    | `purine_count / canonical_base_count`                                                      | DNAKit descriptor_schema_v1                                        |
|  15 | `pyrimidine_count`                         | `composition`           | `count`       | `count(C)+count(T)`                                                                        | DNAKit descriptor_schema_v1                                        |
|  16 | `pyrimidine_fraction`                      | `composition`           | `fraction`    | `pyrimidine_count / canonical_base_count`                                                  | DNAKit descriptor_schema_v1                                        |
|  17 | `amino_count`                              | `composition`           | `count`       | `count(A)+count(C)`                                                                        | DNAKit descriptor_schema_v1                                        |
|  18 | `amino_fraction`                           | `composition`           | `fraction`    | `amino_count / canonical_base_count`                                                       | DNAKit descriptor_schema_v1                                        |
|  19 | `keto_count`                               | `composition`           | `count`       | `count(G)+count(T)`                                                                        | DNAKit descriptor_schema_v1                                        |
|  20 | `keto_fraction`                            | `composition`           | `fraction`    | `keto_count / canonical_base_count`                                                        | DNAKit descriptor_schema_v1                                        |
|  21 | `weak_count`                               | `composition`           | `count`       | `count(A)+count(T)`                                                                        | DNAKit descriptor_schema_v1                                        |
|  22 | `weak_fraction`                            | `composition`           | `fraction`    | `weak_count / canonical_base_count`                                                        | DNAKit descriptor_schema_v1                                        |
|  23 | `strong_count`                             | `composition`           | `count`       | `count(C)+count(G)`                                                                        | DNAKit descriptor_schema_v1                                        |
|  24 | `strong_fraction`                          | `composition`           | `fraction`    | `strong_count / canonical_base_count`                                                      | DNAKit descriptor_schema_v1                                        |
|  25 | `purine_pyrimidine_skew`                   | `composition`           | `ratio`       | `(purine_count-pyrimidine_count)/(purine_count+pyrimidine_count)`                          | DNAKit descriptor_schema_v1                                        |
|  26 | `amino_keto_skew`                          | `composition`           | `ratio`       | `(amino_count-keto_count)/(amino_count+keto_count)`                                        | DNAKit descriptor_schema_v1                                        |
|  27 | `weak_strong_skew`                         | `composition`           | `ratio`       | `(weak_count-strong_count)/(weak_count+strong_count)`                                      | DNAKit descriptor_schema_v1                                        |
|  28 | `gc_at_ratio`                              | `composition`           | `ratio`       | `strong_count / weak_count`                                                                | DNAKit descriptor_schema_v1                                        |
|  29 | `k1_A_frequency`                           | `kmer`                  | `fraction`    | `overlapping count(A) / valid canonical 1-mer positions`                                   | DNAKit exact overlapping k-mer definition                          |
|  30 | `k1_C_frequency`                           | `kmer`                  | `fraction`    | `overlapping count(C) / valid canonical 1-mer positions`                                   | DNAKit exact overlapping k-mer definition                          |
|  31 | `k1_G_frequency`                           | `kmer`                  | `fraction`    | `overlapping count(G) / valid canonical 1-mer positions`                                   | DNAKit exact overlapping k-mer definition                          |
|  32 | `k1_T_frequency`                           | `kmer`                  | `fraction`    | `overlapping count(T) / valid canonical 1-mer positions`                                   | DNAKit exact overlapping k-mer definition                          |
|  33 | `k2_AA_frequency`                          | `kmer`                  | `fraction`    | `overlapping count(AA) / valid canonical 2-mer positions`                                  | DNAKit exact overlapping k-mer definition                          |
|  34 | `k2_AC_frequency`                          | `kmer`                  | `fraction`    | `overlapping count(AC) / valid canonical 2-mer positions`                                  | DNAKit exact overlapping k-mer definition                          |
|  35 | `k2_AG_frequency`                          | `kmer`                  | `fraction`    | `overlapping count(AG) / valid canonical 2-mer positions`                                  | DNAKit exact overlapping k-mer definition                          |
|  36 | `k2_AT_frequency`                          | `kmer`                  | `fraction`    | `overlapping count(AT) / valid canonical 2-mer positions`                                  | DNAKit exact overlapping k-mer definition                          |
|  37 | `k2_CA_frequency`                          | `kmer`                  | `fraction`    | `overlapping count(CA) / valid canonical 2-mer positions`                                  | DNAKit exact overlapping k-mer definition                          |
|  38 | `k2_CC_frequency`                          | `kmer`                  | `fraction`    | `overlapping count(CC) / valid canonical 2-mer positions`                                  | DNAKit exact overlapping k-mer definition                          |
|  39 | `k2_CG_frequency`                          | `kmer`                  | `fraction`    | `overlapping count(CG) / valid canonical 2-mer positions`                                  | DNAKit exact overlapping k-mer definition                          |
|  40 | `k2_CT_frequency`                          | `kmer`                  | `fraction`    | `overlapping count(CT) / valid canonical 2-mer positions`                                  | DNAKit exact overlapping k-mer definition                          |
|  41 | `k2_GA_frequency`                          | `kmer`                  | `fraction`    | `overlapping count(GA) / valid canonical 2-mer positions`                                  | DNAKit exact overlapping k-mer definition                          |
|  42 | `k2_GC_frequency`                          | `kmer`                  | `fraction`    | `overlapping count(GC) / valid canonical 2-mer positions`                                  | DNAKit exact overlapping k-mer definition                          |
|  43 | `k2_GG_frequency`                          | `kmer`                  | `fraction`    | `overlapping count(GG) / valid canonical 2-mer positions`                                  | DNAKit exact overlapping k-mer definition                          |
|  44 | `k2_GT_frequency`                          | `kmer`                  | `fraction`    | `overlapping count(GT) / valid canonical 2-mer positions`                                  | DNAKit exact overlapping k-mer definition                          |
|  45 | `k2_TA_frequency`                          | `kmer`                  | `fraction`    | `overlapping count(TA) / valid canonical 2-mer positions`                                  | DNAKit exact overlapping k-mer definition                          |
|  46 | `k2_TC_frequency`                          | `kmer`                  | `fraction`    | `overlapping count(TC) / valid canonical 2-mer positions`                                  | DNAKit exact overlapping k-mer definition                          |
|  47 | `k2_TG_frequency`                          | `kmer`                  | `fraction`    | `overlapping count(TG) / valid canonical 2-mer positions`                                  | DNAKit exact overlapping k-mer definition                          |
|  48 | `k2_TT_frequency`                          | `kmer`                  | `fraction`    | `overlapping count(TT) / valid canonical 2-mer positions`                                  | DNAKit exact overlapping k-mer definition                          |
|  49 | `k3_AAA_frequency`                         | `kmer`                  | `fraction`    | `overlapping count(AAA) / valid canonical 3-mer positions`                                 | DNAKit exact overlapping k-mer definition                          |
|  50 | `k3_AAC_frequency`                         | `kmer`                  | `fraction`    | `overlapping count(AAC) / valid canonical 3-mer positions`                                 | DNAKit exact overlapping k-mer definition                          |
|  51 | `k3_AAG_frequency`                         | `kmer`                  | `fraction`    | `overlapping count(AAG) / valid canonical 3-mer positions`                                 | DNAKit exact overlapping k-mer definition                          |
|  52 | `k3_AAT_frequency`                         | `kmer`                  | `fraction`    | `overlapping count(AAT) / valid canonical 3-mer positions`                                 | DNAKit exact overlapping k-mer definition                          |
|  53 | `k3_ACA_frequency`                         | `kmer`                  | `fraction`    | `overlapping count(ACA) / valid canonical 3-mer positions`                                 | DNAKit exact overlapping k-mer definition                          |
|  54 | `k3_ACC_frequency`                         | `kmer`                  | `fraction`    | `overlapping count(ACC) / valid canonical 3-mer positions`                                 | DNAKit exact overlapping k-mer definition                          |
|  55 | `k3_ACG_frequency`                         | `kmer`                  | `fraction`    | `overlapping count(ACG) / valid canonical 3-mer positions`                                 | DNAKit exact overlapping k-mer definition                          |
|  56 | `k3_ACT_frequency`                         | `kmer`                  | `fraction`    | `overlapping count(ACT) / valid canonical 3-mer positions`                                 | DNAKit exact overlapping k-mer definition                          |
|  57 | `k3_AGA_frequency`                         | `kmer`                  | `fraction`    | `overlapping count(AGA) / valid canonical 3-mer positions`                                 | DNAKit exact overlapping k-mer definition                          |
|  58 | `k3_AGC_frequency`                         | `kmer`                  | `fraction`    | `overlapping count(AGC) / valid canonical 3-mer positions`                                 | DNAKit exact overlapping k-mer definition                          |
|  59 | `k3_AGG_frequency`                         | `kmer`                  | `fraction`    | `overlapping count(AGG) / valid canonical 3-mer positions`                                 | DNAKit exact overlapping k-mer definition                          |
|  60 | `k3_AGT_frequency`                         | `kmer`                  | `fraction`    | `overlapping count(AGT) / valid canonical 3-mer positions`                                 | DNAKit exact overlapping k-mer definition                          |
|  61 | `k3_ATA_frequency`                         | `kmer`                  | `fraction`    | `overlapping count(ATA) / valid canonical 3-mer positions`                                 | DNAKit exact overlapping k-mer definition                          |
|  62 | `k3_ATC_frequency`                         | `kmer`                  | `fraction`    | `overlapping count(ATC) / valid canonical 3-mer positions`                                 | DNAKit exact overlapping k-mer definition                          |
|  63 | `k3_ATG_frequency`                         | `kmer`                  | `fraction`    | `overlapping count(ATG) / valid canonical 3-mer positions`                                 | DNAKit exact overlapping k-mer definition                          |
|  64 | `k3_ATT_frequency`                         | `kmer`                  | `fraction`    | `overlapping count(ATT) / valid canonical 3-mer positions`                                 | DNAKit exact overlapping k-mer definition                          |
|  65 | `k3_CAA_frequency`                         | `kmer`                  | `fraction`    | `overlapping count(CAA) / valid canonical 3-mer positions`                                 | DNAKit exact overlapping k-mer definition                          |
|  66 | `k3_CAC_frequency`                         | `kmer`                  | `fraction`    | `overlapping count(CAC) / valid canonical 3-mer positions`                                 | DNAKit exact overlapping k-mer definition                          |
|  67 | `k3_CAG_frequency`                         | `kmer`                  | `fraction`    | `overlapping count(CAG) / valid canonical 3-mer positions`                                 | DNAKit exact overlapping k-mer definition                          |
|  68 | `k3_CAT_frequency`                         | `kmer`                  | `fraction`    | `overlapping count(CAT) / valid canonical 3-mer positions`                                 | DNAKit exact overlapping k-mer definition                          |
|  69 | `k3_CCA_frequency`                         | `kmer`                  | `fraction`    | `overlapping count(CCA) / valid canonical 3-mer positions`                                 | DNAKit exact overlapping k-mer definition                          |
|  70 | `k3_CCC_frequency`                         | `kmer`                  | `fraction`    | `overlapping count(CCC) / valid canonical 3-mer positions`                                 | DNAKit exact overlapping k-mer definition                          |
|  71 | `k3_CCG_frequency`                         | `kmer`                  | `fraction`    | `overlapping count(CCG) / valid canonical 3-mer positions`                                 | DNAKit exact overlapping k-mer definition                          |
|  72 | `k3_CCT_frequency`                         | `kmer`                  | `fraction`    | `overlapping count(CCT) / valid canonical 3-mer positions`                                 | DNAKit exact overlapping k-mer definition                          |
|  73 | `k3_CGA_frequency`                         | `kmer`                  | `fraction`    | `overlapping count(CGA) / valid canonical 3-mer positions`                                 | DNAKit exact overlapping k-mer definition                          |
|  74 | `k3_CGC_frequency`                         | `kmer`                  | `fraction`    | `overlapping count(CGC) / valid canonical 3-mer positions`                                 | DNAKit exact overlapping k-mer definition                          |
|  75 | `k3_CGG_frequency`                         | `kmer`                  | `fraction`    | `overlapping count(CGG) / valid canonical 3-mer positions`                                 | DNAKit exact overlapping k-mer definition                          |
|  76 | `k3_CGT_frequency`                         | `kmer`                  | `fraction`    | `overlapping count(CGT) / valid canonical 3-mer positions`                                 | DNAKit exact overlapping k-mer definition                          |
|  77 | `k3_CTA_frequency`                         | `kmer`                  | `fraction`    | `overlapping count(CTA) / valid canonical 3-mer positions`                                 | DNAKit exact overlapping k-mer definition                          |
|  78 | `k3_CTC_frequency`                         | `kmer`                  | `fraction`    | `overlapping count(CTC) / valid canonical 3-mer positions`                                 | DNAKit exact overlapping k-mer definition                          |
|  79 | `k3_CTG_frequency`                         | `kmer`                  | `fraction`    | `overlapping count(CTG) / valid canonical 3-mer positions`                                 | DNAKit exact overlapping k-mer definition                          |
|  80 | `k3_CTT_frequency`                         | `kmer`                  | `fraction`    | `overlapping count(CTT) / valid canonical 3-mer positions`                                 | DNAKit exact overlapping k-mer definition                          |
|  81 | `k3_GAA_frequency`                         | `kmer`                  | `fraction`    | `overlapping count(GAA) / valid canonical 3-mer positions`                                 | DNAKit exact overlapping k-mer definition                          |
|  82 | `k3_GAC_frequency`                         | `kmer`                  | `fraction`    | `overlapping count(GAC) / valid canonical 3-mer positions`                                 | DNAKit exact overlapping k-mer definition                          |
|  83 | `k3_GAG_frequency`                         | `kmer`                  | `fraction`    | `overlapping count(GAG) / valid canonical 3-mer positions`                                 | DNAKit exact overlapping k-mer definition                          |
|  84 | `k3_GAT_frequency`                         | `kmer`                  | `fraction`    | `overlapping count(GAT) / valid canonical 3-mer positions`                                 | DNAKit exact overlapping k-mer definition                          |
|  85 | `k3_GCA_frequency`                         | `kmer`                  | `fraction`    | `overlapping count(GCA) / valid canonical 3-mer positions`                                 | DNAKit exact overlapping k-mer definition                          |
|  86 | `k3_GCC_frequency`                         | `kmer`                  | `fraction`    | `overlapping count(GCC) / valid canonical 3-mer positions`                                 | DNAKit exact overlapping k-mer definition                          |
|  87 | `k3_GCG_frequency`                         | `kmer`                  | `fraction`    | `overlapping count(GCG) / valid canonical 3-mer positions`                                 | DNAKit exact overlapping k-mer definition                          |
|  88 | `k3_GCT_frequency`                         | `kmer`                  | `fraction`    | `overlapping count(GCT) / valid canonical 3-mer positions`                                 | DNAKit exact overlapping k-mer definition                          |
|  89 | `k3_GGA_frequency`                         | `kmer`                  | `fraction`    | `overlapping count(GGA) / valid canonical 3-mer positions`                                 | DNAKit exact overlapping k-mer definition                          |
|  90 | `k3_GGC_frequency`                         | `kmer`                  | `fraction`    | `overlapping count(GGC) / valid canonical 3-mer positions`                                 | DNAKit exact overlapping k-mer definition                          |
|  91 | `k3_GGG_frequency`                         | `kmer`                  | `fraction`    | `overlapping count(GGG) / valid canonical 3-mer positions`                                 | DNAKit exact overlapping k-mer definition                          |
|  92 | `k3_GGT_frequency`                         | `kmer`                  | `fraction`    | `overlapping count(GGT) / valid canonical 3-mer positions`                                 | DNAKit exact overlapping k-mer definition                          |
|  93 | `k3_GTA_frequency`                         | `kmer`                  | `fraction`    | `overlapping count(GTA) / valid canonical 3-mer positions`                                 | DNAKit exact overlapping k-mer definition                          |
|  94 | `k3_GTC_frequency`                         | `kmer`                  | `fraction`    | `overlapping count(GTC) / valid canonical 3-mer positions`                                 | DNAKit exact overlapping k-mer definition                          |
|  95 | `k3_GTG_frequency`                         | `kmer`                  | `fraction`    | `overlapping count(GTG) / valid canonical 3-mer positions`                                 | DNAKit exact overlapping k-mer definition                          |
|  96 | `k3_GTT_frequency`                         | `kmer`                  | `fraction`    | `overlapping count(GTT) / valid canonical 3-mer positions`                                 | DNAKit exact overlapping k-mer definition                          |
|  97 | `k3_TAA_frequency`                         | `kmer`                  | `fraction`    | `overlapping count(TAA) / valid canonical 3-mer positions`                                 | DNAKit exact overlapping k-mer definition                          |
|  98 | `k3_TAC_frequency`                         | `kmer`                  | `fraction`    | `overlapping count(TAC) / valid canonical 3-mer positions`                                 | DNAKit exact overlapping k-mer definition                          |
|  99 | `k3_TAG_frequency`                         | `kmer`                  | `fraction`    | `overlapping count(TAG) / valid canonical 3-mer positions`                                 | DNAKit exact overlapping k-mer definition                          |
| 100 | `k3_TAT_frequency`                         | `kmer`                  | `fraction`    | `overlapping count(TAT) / valid canonical 3-mer positions`                                 | DNAKit exact overlapping k-mer definition                          |
| 101 | `k3_TCA_frequency`                         | `kmer`                  | `fraction`    | `overlapping count(TCA) / valid canonical 3-mer positions`                                 | DNAKit exact overlapping k-mer definition                          |
| 102 | `k3_TCC_frequency`                         | `kmer`                  | `fraction`    | `overlapping count(TCC) / valid canonical 3-mer positions`                                 | DNAKit exact overlapping k-mer definition                          |
| 103 | `k3_TCG_frequency`                         | `kmer`                  | `fraction`    | `overlapping count(TCG) / valid canonical 3-mer positions`                                 | DNAKit exact overlapping k-mer definition                          |
| 104 | `k3_TCT_frequency`                         | `kmer`                  | `fraction`    | `overlapping count(TCT) / valid canonical 3-mer positions`                                 | DNAKit exact overlapping k-mer definition                          |
| 105 | `k3_TGA_frequency`                         | `kmer`                  | `fraction`    | `overlapping count(TGA) / valid canonical 3-mer positions`                                 | DNAKit exact overlapping k-mer definition                          |
| 106 | `k3_TGC_frequency`                         | `kmer`                  | `fraction`    | `overlapping count(TGC) / valid canonical 3-mer positions`                                 | DNAKit exact overlapping k-mer definition                          |
| 107 | `k3_TGG_frequency`                         | `kmer`                  | `fraction`    | `overlapping count(TGG) / valid canonical 3-mer positions`                                 | DNAKit exact overlapping k-mer definition                          |
| 108 | `k3_TGT_frequency`                         | `kmer`                  | `fraction`    | `overlapping count(TGT) / valid canonical 3-mer positions`                                 | DNAKit exact overlapping k-mer definition                          |
| 109 | `k3_TTA_frequency`                         | `kmer`                  | `fraction`    | `overlapping count(TTA) / valid canonical 3-mer positions`                                 | DNAKit exact overlapping k-mer definition                          |
| 110 | `k3_TTC_frequency`                         | `kmer`                  | `fraction`    | `overlapping count(TTC) / valid canonical 3-mer positions`                                 | DNAKit exact overlapping k-mer definition                          |
| 111 | `k3_TTG_frequency`                         | `kmer`                  | `fraction`    | `overlapping count(TTG) / valid canonical 3-mer positions`                                 | DNAKit exact overlapping k-mer definition                          |
| 112 | `k3_TTT_frequency`                         | `kmer`                  | `fraction`    | `overlapping count(TTT) / valid canonical 3-mer positions`                                 | DNAKit exact overlapping k-mer definition                          |
| 113 | `gc_skew`                                  | `skew_cpg`              | `ratio`       | `(count(G)-count(C))/(count(G)+count(C))`                                                  | DNAKit descriptor_schema_v1                                        |
| 114 | `at_skew`                                  | `skew_cpg`              | `ratio`       | `(count(A)-count(T))/(count(A)+count(T))`                                                  | DNAKit descriptor_schema_v1                                        |
| 115 | `cpg_count`                                | `skew_cpg`              | `count`       | `overlapping count(CG)`                                                                    | DNAKit descriptor_schema_v1                                        |
| 116 | `cpg_density`                              | `skew_cpg`              | `fraction`    | `count(CG) / valid canonical dinucleotide positions`                                       | DNAKit descriptor_schema_v1                                        |
| 117 | `cpg_observed_expected`                    | `skew_cpg`              | `ratio`       | `count(CG)*canonical_base_count/(count(C)*count(G))`                                       | DNAKit descriptor_schema_v1                                        |
| 118 | `gpc_count`                                | `skew_cpg`              | `count`       | `overlapping count(GC)`                                                                    | DNAKit descriptor_schema_v1                                        |
| 119 | `gpc_density`                              | `skew_cpg`              | `fraction`    | `count(GC) / valid canonical dinucleotide positions`                                       | DNAKit descriptor_schema_v1                                        |
| 120 | `cpg_gpc_ratio`                            | `skew_cpg`              | `ratio`       | `count(CG) / count(GC)`                                                                    | DNAKit descriptor_schema_v1                                        |
| 121 | `cumulative_gc_skew_max`                   | `skew_cpg`              | `count`       | `max prefix cumulative score where G=+1,C=-1,A/T=0`                                        | DNAKit descriptor_schema_v1                                        |
| 122 | `cumulative_gc_skew_min`                   | `skew_cpg`              | `count`       | `min prefix cumulative score where G=+1,C=-1,A/T=0`                                        | DNAKit descriptor_schema_v1                                        |
| 123 | `cumulative_gc_skew_range`                 | `skew_cpg`              | `count`       | `cumulative_gc_skew_max - cumulative_gc_skew_min`                                          | DNAKit descriptor_schema_v1                                        |
| 124 | `cumulative_at_skew_max`                   | `skew_cpg`              | `count`       | `max prefix cumulative score where A=+1,T=-1,C/G=0`                                        | DNAKit descriptor_schema_v1                                        |
| 125 | `cumulative_at_skew_min`                   | `skew_cpg`              | `count`       | `min prefix cumulative score where A=+1,T=-1,C/G=0`                                        | DNAKit descriptor_schema_v1                                        |
| 126 | `cumulative_at_skew_range`                 | `skew_cpg`              | `count`       | `cumulative_at_skew_max - cumulative_at_skew_min`                                          | DNAKit descriptor_schema_v1                                        |
| 127 | `dinucleotide_rc_total_variation`          | `skew_cpg`              | `fraction`    | `0.5*sum_xy(abs(f_xy-f_reverse_complement(xy)))`                                           | DNAKit descriptor_schema_v1                                        |
| 128 | `mono_chargaff_l1_distance`                | `skew_cpg`              | `fraction`    | `abs(f_A-f_T)+abs(f_C-f_G)`                                                                | DNAKit descriptor_schema_v1                                        |
| 129 | `shannon_entropy_k1_bits`                  | `complexity`            | `bits`        | `-sum(p(1-mer)*log2(p(1-mer)))`                                                            | Shannon 1948; DOI 10.1002/j.1538-7305.1948.tb01338.x               |
| 130 | `shannon_entropy_k2_bits`                  | `complexity`            | `bits`        | `-sum(p(2-mer)*log2(p(2-mer)))`                                                            | Shannon 1948; DOI 10.1002/j.1538-7305.1948.tb01338.x               |
| 131 | `shannon_entropy_k3_bits`                  | `complexity`            | `bits`        | `-sum(p(3-mer)*log2(p(3-mer)))`                                                            | Shannon 1948; DOI 10.1002/j.1538-7305.1948.tb01338.x               |
| 132 | `normalized_entropy_k1`                    | `complexity`            | `fraction`    | `shannon_entropy_k1_bits / log2(4**1)`                                                     | Shannon 1948; DOI 10.1002/j.1538-7305.1948.tb01338.x               |
| 133 | `normalized_entropy_k2`                    | `complexity`            | `fraction`    | `shannon_entropy_k2_bits / log2(4**2)`                                                     | Shannon 1948; DOI 10.1002/j.1538-7305.1948.tb01338.x               |
| 134 | `normalized_entropy_k3`                    | `complexity`            | `fraction`    | `shannon_entropy_k3_bits / log2(4**3)`                                                     | Shannon 1948; DOI 10.1002/j.1538-7305.1948.tb01338.x               |
| 135 | `linguistic_complexity_k2`                 | `complexity`            | `fraction`    | `unique 2-mers / min(4**2, valid 2-mer positions)`                                         | Observed/possible k-word linguistic complexity                     |
| 136 | `linguistic_complexity_k3`                 | `complexity`            | `fraction`    | `unique 3-mers / min(4**3, valid 3-mer positions)`                                         | Observed/possible k-word linguistic complexity                     |
| 137 | `linguistic_complexity_k4`                 | `complexity`            | `fraction`    | `unique 4-mers / min(4**4, valid 4-mer positions)`                                         | Observed/possible k-word linguistic complexity                     |
| 138 | `linguistic_complexity_k5`                 | `complexity`            | `fraction`    | `unique 5-mers / min(4**5, valid 5-mer positions)`                                         | Observed/possible k-word linguistic complexity                     |
| 139 | `linguistic_complexity_k6`                 | `complexity`            | `fraction`    | `unique 6-mers / min(4**6, valid 6-mer positions)`                                         | Observed/possible k-word linguistic complexity                     |
| 140 | `linguistic_complexity_product_k1_k6`      | `complexity`            | `fraction`    | `product of defined linguistic_complexity_k values for k=1..6`                             | Observed/possible k-word linguistic complexity                     |
| 141 | `lz76_complexity`                          | `complexity`            | `count`       | `number of phrases in exhaustive Lempel-Ziv 1976 parsing`                                  | Lempel and Ziv 1976; DOI 10.1109/TIT.1976.1055501                  |
| 142 | `normalized_lz76_complexity`               | `complexity`            | `ratio`       | `lz76_complexity*log_base4(canonical_base_count)/canonical_base_count`                     | Lempel and Ziv 1976; DOI 10.1109/TIT.1976.1055501                  |
| 143 | `longest_homopolymer_nt`                   | `complexity`            | `nt`          | `max canonical homopolymer run`                                                            | DNAKit descriptor_schema_v1                                        |
| 144 | `longest_homopolymer_a_nt`                 | `complexity`            | `nt`          | `max homopolymer run of A`                                                                 | DNAKit descriptor_schema_v1                                        |
| 145 | `longest_homopolymer_c_nt`                 | `complexity`            | `nt`          | `max homopolymer run of C`                                                                 | DNAKit descriptor_schema_v1                                        |
| 146 | `longest_homopolymer_g_nt`                 | `complexity`            | `nt`          | `max homopolymer run of G`                                                                 | DNAKit descriptor_schema_v1                                        |
| 147 | `longest_homopolymer_t_nt`                 | `complexity`            | `nt`          | `max homopolymer run of T`                                                                 | DNAKit descriptor_schema_v1                                        |
| 148 | `exact_tandem_repeat_coverage_fraction`    | `complexity`            | `fraction`    | `union bases covered by exact tandem repeats / canonical_base_count`                       | DNAKit exact tandem repeat scanner; units 1..20; minimum repeats 2 |
| 149 | `frame0_codon_count`                       | `coding`                | `count`       | `number of valid forward frame-0 codons`                                                   | NCBI standard genetic code table 1; DNAKit ORF rules               |
| 150 | `frame0_unique_codon_count`                | `coding`                | `count`       | `number of distinct forward frame-0 codons`                                                | NCBI standard genetic code table 1; DNAKit ORF rules               |
| 151 | `frame0_start_codon_count`                 | `coding`                | `count`       | `count(ATG) in forward frame 0`                                                            | NCBI standard genetic code table 1; DNAKit ORF rules               |
| 152 | `frame0_stop_codon_count`                  | `coding`                | `count`       | `count(TAA,TAG,TGA) in forward frame 0`                                                    | NCBI standard genetic code table 1; DNAKit ORF rules               |
| 153 | `frame0_start_codon_fraction`              | `coding`                | `fraction`    | `frame0_start_codon_count / frame0_codon_count`                                            | NCBI standard genetic code table 1; DNAKit ORF rules               |
| 154 | `frame0_stop_codon_fraction`               | `coding`                | `fraction`    | `frame0_stop_codon_count / frame0_codon_count`                                             | NCBI standard genetic code table 1; DNAKit ORF rules               |
| 155 | `frame0_codon_entropy_bits`                | `coding`                | `bits`        | `-sum(frame0 codon frequency*log2(frequency))`                                             | NCBI standard genetic code table 1; DNAKit ORF rules               |
| 156 | `frame0_effective_number_of_codons`        | `coding`                | `count`       | `2**frame0_codon_entropy_bits`                                                             | NCBI standard genetic code table 1; DNAKit ORF rules               |
| 157 | `frame0_gc1_fraction`                      | `coding`                | `fraction`    | `GC bases at position 1 / valid frame-0 codons`                                            | NCBI standard genetic code table 1; DNAKit ORF rules               |
| 158 | `frame0_gc2_fraction`                      | `coding`                | `fraction`    | `GC bases at position 2 / valid frame-0 codons`                                            | NCBI standard genetic code table 1; DNAKit ORF rules               |
| 159 | `frame0_gc3_fraction`                      | `coding`                | `fraction`    | `GC bases at position 3 / valid frame-0 codons`                                            | NCBI standard genetic code table 1; DNAKit ORF rules               |
| 160 | `six_frame_complete_orf_count`             | `coding`                | `count`       | `complete start-to-next-stop ORFs across three frames on both strands`                     | NCBI standard genetic code table 1; DNAKit ORF rules               |
| 161 | `six_frame_forward_complete_orf_count`     | `coding`                | `count`       | `complete ORFs across three forward frames`                                                | NCBI standard genetic code table 1; DNAKit ORF rules               |
| 162 | `six_frame_reverse_complete_orf_count`     | `coding`                | `count`       | `complete ORFs across three reverse-complement frames`                                     | NCBI standard genetic code table 1; DNAKit ORF rules               |
| 163 | `six_frame_longest_complete_orf_nt`        | `coding`                | `nt`          | `maximum complete six-frame ORF length including terminal stop`                            | NCBI standard genetic code table 1; DNAKit ORF rules               |
| 164 | `six_frame_complete_orf_coverage_fraction` | `coding`                | `fraction`    | `union of symbol positions covered by complete six-frame ORFs / canonical_base_count`      | NCBI standard genetic code table 1; DNAKit ORF rules               |
| 165 | `mw_ss_oh_da`                              | `physicochemical`       | `Da`          | `anhydrous mass of one ssDNA strand with 5-prime OH`                                       | DNAKit anhydrous DNA residue mass table v1                         |
| 166 | `mw_ss_5p_phosphate_da`                    | `physicochemical`       | `Da`          | `anhydrous mass of one ssDNA strand with 5-prime phosphate`                                | DNAKit anhydrous DNA residue mass table v1                         |
| 167 | `mw_ds_oh_da`                              | `physicochemical`       | `Da`          | `anhydrous mass of sequence plus complete reverse complement with 5-prime OH`              | DNAKit anhydrous DNA residue mass table v1                         |
| 168 | `mw_ds_5p_phosphate_da`                    | `physicochemical`       | `Da`          | `anhydrous mass of sequence plus complete reverse complement; both 5-prime phosphorylated` | DNAKit anhydrous DNA residue mass table v1                         |
| 169 | `epsilon260_ss_m_inverse_cm_inverse`       | `physicochemical`       | `M^-1 cm^-1`  | `nearest-neighbor epsilon260 pair sum minus internal-base sum`                             | Warshaw-Tinoco 1966 and Cantor-Warshaw-Shapiro 1970                |
| 170 | `nmol_per_a260_1ml_1cm`                    | `physicochemical`       | `nmol`        | `1e6 / epsilon260 for A260=1, volume=1 mL, path=1 cm`                                      | Warshaw-Tinoco 1966 and Cantor-Warshaw-Shapiro 1970                |
| 171 | `ug_per_a260_1ml_1cm`                      | `physicochemical`       | `ug`          | `1000*mw_ss_oh_da/epsilon260 for A260=1, volume=1 mL, path=1 cm`                           | Warshaw-Tinoco 1966 and Cantor-Warshaw-Shapiro 1970                |
| 172 | `tm_wallace_c`                             | `physicochemical`       | `degree C`    | `2*(A+T)+4*(G+C)`                                                                          | Wallace short-oligo 2AT+4GC rule                                   |
| 173 | `stacking_delta_h_kcal_per_mol`            | `physicochemical`       | `kcal/mol`    | `sum SantaLucia nearest-neighbor stacking delta H`                                         | SantaLucia 1998; DOI 10.1073/pnas.95.4.1460                        |
| 174 | `stacking_delta_s_cal_per_k_mol_k`         | `physicochemical`       | `cal/(K mol)` | `sum SantaLucia nearest-neighbor stacking delta S`                                         | SantaLucia 1998; DOI 10.1073/pnas.95.4.1460                        |
| 175 | `stacking_delta_g37_kcal_per_mol`          | `physicochemical`       | `kcal/mol`    | `stacking delta H - 310.15*stacking delta S/1000`                                          | SantaLucia 1998; DOI 10.1073/pnas.95.4.1460                        |
| 176 | `nn_delta_h_kcal_per_mol`                  | `physicochemical`       | `kcal/mol`    | `SantaLucia complete-duplex delta H with initiation and symmetry`                          | SantaLucia 1998; DOI 10.1073/pnas.95.4.1460                        |
| 177 | `nn_delta_s_cal_per_mol_k`                 | `physicochemical`       | `cal/(K mol)` | `SantaLucia complete-duplex delta S with initiation, symmetry, and salt`                   | SantaLucia 1998; DOI 10.1073/pnas.95.4.1460                        |
| 178 | `nn_delta_g37_kcal_per_mol`                | `physicochemical`       | `kcal/mol`    | `complete-duplex delta H - 310.15*delta S/1000`                                            | SantaLucia 1998; DOI 10.1073/pnas.95.4.1460                        |
| 179 | `nn_tm_c`                                  | `physicochemical`       | `degree C`    | `SantaLucia concentration- and sodium-adjusted Tm`                                         | SantaLucia 1998; DOI 10.1073/pnas.95.4.1460                        |
| 180 | `self_complementary`                       | `physicochemical`       | `boolean`     | `sequence == reverse_complement(sequence)`                                                 | SantaLucia 1998; DOI 10.1073/pnas.95.4.1460                        |
| 181 | `diprodb_twist_mean`                       | `dinucleotide_property` | `degree`      | `population mean of Twist values over valid overlapping dinucleotides`                     | Caller-supplied table; DNAKit bundles no numerical values          |
| 182 | `diprodb_twist_sd`                         | `dinucleotide_property` | `degree`      | `population sd of Twist values over valid overlapping dinucleotides`                       | Caller-supplied table; DNAKit bundles no numerical values          |
| 183 | `diprodb_twist_min`                        | `dinucleotide_property` | `degree`      | `population min of Twist values over valid overlapping dinucleotides`                      | Caller-supplied table; DNAKit bundles no numerical values          |
| 184 | `diprodb_twist_max`                        | `dinucleotide_property` | `degree`      | `population max of Twist values over valid overlapping dinucleotides`                      | Caller-supplied table; DNAKit bundles no numerical values          |
| 185 | `diprodb_tilt_mean`                        | `dinucleotide_property` | `degree`      | `population mean of Tilt values over valid overlapping dinucleotides`                      | Caller-supplied table; DNAKit bundles no numerical values          |
| 186 | `diprodb_tilt_sd`                          | `dinucleotide_property` | `degree`      | `population sd of Tilt values over valid overlapping dinucleotides`                        | Caller-supplied table; DNAKit bundles no numerical values          |
| 187 | `diprodb_tilt_min`                         | `dinucleotide_property` | `degree`      | `population min of Tilt values over valid overlapping dinucleotides`                       | Caller-supplied table; DNAKit bundles no numerical values          |
| 188 | `diprodb_tilt_max`                         | `dinucleotide_property` | `degree`      | `population max of Tilt values over valid overlapping dinucleotides`                       | Caller-supplied table; DNAKit bundles no numerical values          |
| 189 | `diprodb_roll_mean`                        | `dinucleotide_property` | `degree`      | `population mean of Roll values over valid overlapping dinucleotides`                      | Caller-supplied table; DNAKit bundles no numerical values          |
| 190 | `diprodb_roll_sd`                          | `dinucleotide_property` | `degree`      | `population sd of Roll values over valid overlapping dinucleotides`                        | Caller-supplied table; DNAKit bundles no numerical values          |
| 191 | `diprodb_roll_min`                         | `dinucleotide_property` | `degree`      | `population min of Roll values over valid overlapping dinucleotides`                       | Caller-supplied table; DNAKit bundles no numerical values          |
| 192 | `diprodb_roll_max`                         | `dinucleotide_property` | `degree`      | `population max of Roll values over valid overlapping dinucleotides`                       | Caller-supplied table; DNAKit bundles no numerical values          |
| 193 | `diprodb_shift_mean`                       | `dinucleotide_property` | `angstrom`    | `population mean of Shift values over valid overlapping dinucleotides`                     | Caller-supplied table; DNAKit bundles no numerical values          |
| 194 | `diprodb_shift_sd`                         | `dinucleotide_property` | `angstrom`    | `population sd of Shift values over valid overlapping dinucleotides`                       | Caller-supplied table; DNAKit bundles no numerical values          |
| 195 | `diprodb_shift_min`                        | `dinucleotide_property` | `angstrom`    | `population min of Shift values over valid overlapping dinucleotides`                      | Caller-supplied table; DNAKit bundles no numerical values          |
| 196 | `diprodb_shift_max`                        | `dinucleotide_property` | `angstrom`    | `population max of Shift values over valid overlapping dinucleotides`                      | Caller-supplied table; DNAKit bundles no numerical values          |
| 197 | `diprodb_slide_mean`                       | `dinucleotide_property` | `angstrom`    | `population mean of Slide values over valid overlapping dinucleotides`                     | Caller-supplied table; DNAKit bundles no numerical values          |
| 198 | `diprodb_slide_sd`                         | `dinucleotide_property` | `angstrom`    | `population sd of Slide values over valid overlapping dinucleotides`                       | Caller-supplied table; DNAKit bundles no numerical values          |
| 199 | `diprodb_slide_min`                        | `dinucleotide_property` | `angstrom`    | `population min of Slide values over valid overlapping dinucleotides`                      | Caller-supplied table; DNAKit bundles no numerical values          |
| 200 | `diprodb_slide_max`                        | `dinucleotide_property` | `angstrom`    | `population max of Slide values over valid overlapping dinucleotides`                      | Caller-supplied table; DNAKit bundles no numerical values          |
| 201 | `diprodb_rise_mean`                        | `dinucleotide_property` | `angstrom`    | `population mean of Rise values over valid overlapping dinucleotides`                      | Caller-supplied table; DNAKit bundles no numerical values          |
| 202 | `diprodb_rise_sd`                          | `dinucleotide_property` | `angstrom`    | `population sd of Rise values over valid overlapping dinucleotides`                        | Caller-supplied table; DNAKit bundles no numerical values          |
| 203 | `diprodb_rise_min`                         | `dinucleotide_property` | `angstrom`    | `population min of Rise values over valid overlapping dinucleotides`                       | Caller-supplied table; DNAKit bundles no numerical values          |
| 204 | `diprodb_rise_max`                         | `dinucleotide_property` | `angstrom`    | `population max of Rise values over valid overlapping dinucleotides`                       | Caller-supplied table; DNAKit bundles no numerical values          |
| 205 | `diprodb_bend_mean`                        | `dinucleotide_property` | `degree`      | `population mean of Bend values over valid overlapping dinucleotides`                      | Caller-supplied table; DNAKit bundles no numerical values          |
| 206 | `diprodb_bend_sd`                          | `dinucleotide_property` | `degree`      | `population sd of Bend values over valid overlapping dinucleotides`                        | Caller-supplied table; DNAKit bundles no numerical values          |
| 207 | `diprodb_bend_min`                         | `dinucleotide_property` | `degree`      | `population min of Bend values over valid overlapping dinucleotides`                       | Caller-supplied table; DNAKit bundles no numerical values          |
| 208 | `diprodb_bend_max`                         | `dinucleotide_property` | `degree`      | `population max of Bend values over valid overlapping dinucleotides`                       | Caller-supplied table; DNAKit bundles no numerical values          |
| 209 | `diprodb_inclination_mean`                 | `dinucleotide_property` | `degree`      | `population mean of Inclination values over valid overlapping dinucleotides`               | Caller-supplied table; DNAKit bundles no numerical values          |
| 210 | `diprodb_inclination_sd`                   | `dinucleotide_property` | `degree`      | `population sd of Inclination values over valid overlapping dinucleotides`                 | Caller-supplied table; DNAKit bundles no numerical values          |
| 211 | `diprodb_inclination_min`                  | `dinucleotide_property` | `degree`      | `population min of Inclination values over valid overlapping dinucleotides`                | Caller-supplied table; DNAKit bundles no numerical values          |
| 212 | `diprodb_inclination_max`                  | `dinucleotide_property` | `degree`      | `population max of Inclination values over valid overlapping dinucleotides`                | Caller-supplied table; DNAKit bundles no numerical values          |
| 213 | `diprodb_direction_mean`                   | `dinucleotide_property` | `degree`      | `population mean of Direction values over valid overlapping dinucleotides`                 | Caller-supplied table; DNAKit bundles no numerical values          |
| 214 | `diprodb_direction_sd`                     | `dinucleotide_property` | `degree`      | `population sd of Direction values over valid overlapping dinucleotides`                   | Caller-supplied table; DNAKit bundles no numerical values          |
| 215 | `diprodb_direction_min`                    | `dinucleotide_property` | `degree`      | `population min of Direction values over valid overlapping dinucleotides`                  | Caller-supplied table; DNAKit bundles no numerical values          |
| 216 | `diprodb_direction_max`                    | `dinucleotide_property` | `degree`      | `population max of Direction values over valid overlapping dinucleotides`                  | Caller-supplied table; DNAKit bundles no numerical values          |
| 217 | `diprodb_propeller_twist_mean`             | `dinucleotide_property` | `degree`      | `population mean of Propeller twist values over valid overlapping dinucleotides`           | Caller-supplied table; DNAKit bundles no numerical values          |
| 218 | `diprodb_propeller_twist_sd`               | `dinucleotide_property` | `degree`      | `population sd of Propeller twist values over valid overlapping dinucleotides`             | Caller-supplied table; DNAKit bundles no numerical values          |
| 219 | `diprodb_propeller_twist_min`              | `dinucleotide_property` | `degree`      | `population min of Propeller twist values over valid overlapping dinucleotides`            | Caller-supplied table; DNAKit bundles no numerical values          |
| 220 | `diprodb_propeller_twist_max`              | `dinucleotide_property` | `degree`      | `population max of Propeller twist values over valid overlapping dinucleotides`            | Caller-supplied table; DNAKit bundles no numerical values          |
| 221 | `diprodb_major_groove_width_mean`          | `dinucleotide_property` | `angstrom`    | `population mean of Major groove width values over valid overlapping dinucleotides`        | Caller-supplied table; DNAKit bundles no numerical values          |
| 222 | `diprodb_major_groove_width_sd`            | `dinucleotide_property` | `angstrom`    | `population sd of Major groove width values over valid overlapping dinucleotides`          | Caller-supplied table; DNAKit bundles no numerical values          |
| 223 | `diprodb_major_groove_width_min`           | `dinucleotide_property` | `angstrom`    | `population min of Major groove width values over valid overlapping dinucleotides`         | Caller-supplied table; DNAKit bundles no numerical values          |
| 224 | `diprodb_major_groove_width_max`           | `dinucleotide_property` | `angstrom`    | `population max of Major groove width values over valid overlapping dinucleotides`         | Caller-supplied table; DNAKit bundles no numerical values          |
| 225 | `diprodb_minor_groove_width_mean`          | `dinucleotide_property` | `angstrom`    | `population mean of Minor groove width values over valid overlapping dinucleotides`        | Caller-supplied table; DNAKit bundles no numerical values          |
| 226 | `diprodb_minor_groove_width_sd`            | `dinucleotide_property` | `angstrom`    | `population sd of Minor groove width values over valid overlapping dinucleotides`          | Caller-supplied table; DNAKit bundles no numerical values          |
| 227 | `diprodb_minor_groove_width_min`           | `dinucleotide_property` | `angstrom`    | `population min of Minor groove width values over valid overlapping dinucleotides`         | Caller-supplied table; DNAKit bundles no numerical values          |
| 228 | `diprodb_minor_groove_width_max`           | `dinucleotide_property` | `angstrom`    | `population max of Minor groove width values over valid overlapping dinucleotides`         | Caller-supplied table; DNAKit bundles no numerical values          |
| 229 | `diprodb_persistence_length_mean`          | `dinucleotide_property` | `nanometer`   | `population mean of Persistence length values over valid overlapping dinucleotides`        | Caller-supplied table; DNAKit bundles no numerical values          |
| 230 | `diprodb_persistence_length_sd`            | `dinucleotide_property` | `nanometer`   | `population sd of Persistence length values over valid overlapping dinucleotides`          | Caller-supplied table; DNAKit bundles no numerical values          |
| 231 | `diprodb_persistence_length_min`           | `dinucleotide_property` | `nanometer`   | `population min of Persistence length values over valid overlapping dinucleotides`         | Caller-supplied table; DNAKit bundles no numerical values          |
| 232 | `diprodb_persistence_length_max`           | `dinucleotide_property` | `nanometer`   | `population max of Persistence length values over valid overlapping dinucleotides`         | Caller-supplied table; DNAKit bundles no numerical values          |
| 233 | `diprodb_stacking_energy_mean`             | `dinucleotide_property` | `kcal/mol`    | `population mean of Stacking energy values over valid overlapping dinucleotides`           | Caller-supplied table; DNAKit bundles no numerical values          |
| 234 | `diprodb_stacking_energy_sd`               | `dinucleotide_property` | `kcal/mol`    | `population sd of Stacking energy values over valid overlapping dinucleotides`             | Caller-supplied table; DNAKit bundles no numerical values          |
| 235 | `diprodb_stacking_energy_min`              | `dinucleotide_property` | `kcal/mol`    | `population min of Stacking energy values over valid overlapping dinucleotides`            | Caller-supplied table; DNAKit bundles no numerical values          |
| 236 | `diprodb_stacking_energy_max` | `dinucleotide_property` | `kcal/mol` | `population max of Stacking energy values over valid overlapping dinucleotides` | Caller-supplied table; DNAKit bundles no numerical values |
| 237 | `diprodb_free_energy_mean` | `dinucleotide_property` | `kcal/mol` | `population mean of Free energy values over valid overlapping dinucleotides` | Caller-supplied table; DNAKit bundles no numerical values |
| 238 | `diprodb_free_energy_sd` | `dinucleotide_property` | `kcal/mol` | `population sd of Free energy values over valid overlapping dinucleotides` | Caller-supplied table; DNAKit bundles no numerical values |
| 239 | `diprodb_free_energy_min` | `dinucleotide_property` | `kcal/mol` | `population min of Free energy values over valid overlapping dinucleotides` | Caller-supplied table; DNAKit bundles no numerical values |
| 240 | `diprodb_free_energy_max` | `dinucleotide_property` | `kcal/mol` | `population max of Free energy values over valid overlapping dinucleotides` | Caller-supplied table; DNAKit bundles no numerical values |

Methods, papers, databases and URLs have been moved to [Acknowledgments and Primary Sources](../../acknowledgements.md#methods-and-references); permissions and user responsibilities can be found in [Third Party Statements](../../acknowledgements.md#third-party-notices).



<span id="all-descriptors"></span>
