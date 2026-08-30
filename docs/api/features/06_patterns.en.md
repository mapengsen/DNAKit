# Sequence function search

Search DNA sequences for patterns such as motifs/PWMs, functional sites, palindromes, inverted repeats, tandem repeats, and microsatellites.

## 1) PAT-001 / SIM-003 motif search

- **Meaning:** A motif is a short sequence pattern with certain characteristics in DNA, which can be understood as a "keyword" in the DNA sequence.
- **Function:** Scan the specified motif according to exact, IUPAC or restricted regular rules, and return the coordinates, matching text and chain direction of each hit, which can be used to locate known sequence patterns.
- **API:** `dnakit.patterns.scan_motif(value[required], motif[required], mode[optional], name[optional], strand[optional], overlapping[optional], merge_strands[optional], max_matches[optional], max_scan_length[optional], max_pattern_length[optional], max_scan_cells[optional])`, `dnakit.patterns.scan_pwm(value[required], pwm[required], threshold[required], background[optional], pseudocount[optional], strand[optional], max_matches[optional], max_scan_length[optional], max_pwm_length[optional], max_score_cells[optional])`, `dnakit.patterns.PWM(name[required], matrix[required])`
- **Input:** Required sequence and motif/regex/PWM; optional pattern, strand, overlap, background, threshold and hit cap.
- **Sample code:**

```python
from dnakit import DNASequence
from dnakit.patterns import scan_motif

result = scan_motif(DNASequence("AAAA"), "AA", strand="forward")
print([(hit.symbol_location.start, hit.symbol_location.end) for hit in result.hits])
```

- **Example results:**

```text
[(0, 2), (1, 3), (2, 4)]
```

## 2) PAT-003 start and stop codons

- **Meaning:** A start codon indicates where protein translation may begin, and a stop codon indicates where protein translation ends.
- **Calculation rules:** Deterministic rules ([FAQ detailed explanation](../../faq.md#pattern-matching-strategy)).
- **Function:** Search the start and stop codons in the six reading frames of both forward and reverse strands, return the type, reading frame and position, and provide candidate sites for ORF and coding region analysis.
- **API:** `dnakit.patterns.scan_codon_sites(value[required], genetic_code[optional], start_codons[optional], stop_codons[optional], strand[optional], max_matches[optional], max_codon_checks[optional])`
- **Input:** Required sequence; optional genetic code table, custom start/stop set, and strand.
- **Sample code:**

```python
from dnakit import DNASequence
from dnakit.patterns import scan_codon_sites

result = scan_codon_sites(DNASequence("ATGAAATAA"), strand="forward")
print([(hit.kind, hit.codon, hit.frame) for hit in result.hits])
```

- **Example results:**

```text
[('start', 'ATG', 1), ('stop', 'TAA', 1), ('stop', 'TGA', 2)]
```

## 3) PAT-004 promoter motif

- **Meaning:** A promoter is a regulatory region that helps genes start to be transcribed; the promoter motif is a common short sequence feature and can be understood as a "start mark" for gene transcription.
- **Calculation rules:** Deterministic rules ([FAQ detailed explanation](../../faq.md#pattern-matching-strategy)).
- **Function:** Scan the built-in consensus sequence or user-provided promoter motif and return candidate positions and chain directions for regular promoter region screening, but does not predict promoter activity.
- **API:** `dnakit.patterns.scan_promoter_motifs(value[required], motifs[optional], strand[optional], max_matches[optional], max_scan_length[optional], max_scan_cells[optional], max_motifs[optional])`
- **Input:** Required sequence; optional motif mapping, strand and resource cap.
- **Sample code:**

```python
from dnakit import DNASequence
from dnakit.patterns import scan_promoter_motifs

result = scan_promoter_motifs(DNASequence("GGTATAATCC"), strand="forward")
print([hit.motif_name for hit in result.hits])
```

- **Example results:**

```text
['bacterial_minus_10_consensus']
```

## 4) PAT-005 TF motif

- **Meaning:** TF motif is a DNA sequence pattern that transcription factors preferentially recognize and bind to, and may be involved in regulating whether genes are expressed and how strongly they are expressed.
- **Calculation rules:** Deterministic rules ([FAQ detailed explanation](../../faq.md#pattern-matching-strategy)).
- **Function:** Use the position weight matrix provided by the user for window-by-window scoring, and return candidate transcription factor binding sites and scores that exceed the threshold, which are used for motif candidate screening rather than actual binding strength prediction.
- **API:** `dnakit.patterns.PWM(name[required], matrix[required])`, `dnakit.patterns.scan_tf_pwm(value[required], tf_name[required], pwm[required], threshold[required], background[optional], pseudocount[optional], strand[optional], max_matches[optional], max_scan_length[optional], max_pwm_length[optional], max_score_cells[optional])`
- **Input:** Required sequence, TF name, PWM and threshold; optional background, pseudocount and strand.
- **Sample code:**

```python
from dnakit import DNASequence
from dnakit.patterns import PWM, scan_tf_pwm

pwm = PWM("input", {"A": [4], "C": [0], "G": [0], "T": [0]})
result = scan_tf_pwm(DNASequence("AA"), "TF-X", pwm, threshold=1.0, strand="forward")
print([hit.motif_name for hit in result.hits])
```

- **Example results:**

```text
['TF-X', 'TF-X']
```

## 5) PAT-006 restriction enzyme site

- **Meaning:** Restriction enzyme sites are specific DNA sequences that restriction enzymes can recognize and cut. They can be understood as the cutting positions of "molecular scissors".
- **Calculation rules:** Deterministic rules ([FAQ detailed explanation](../../faq.md#pattern-matching-strategy)).
- **Function:** Search for recognition sequences based on built-in or custom restriction enzyme definitions, and return forward and reverse strand hits and cutting coordinates, which are used to plan enzyme digestion experiments and determine fragment boundaries.
- **API:** `dnakit.patterns.scan_restriction_sites(value[required], enzymes[required], max_matches[optional], max_scan_length[optional], max_scan_cells[optional], max_enzymes[optional])`, `dnakit.patterns.RestrictionEnzyme(name[required], recognition_sequence[required], top_cut_offset[required], bottom_cut_offset[required], source[optional])`
- **Input:** Required sequence and enzyme name/definition list; optional scan and hit cap.
- **Sample code:**

```python
from dnakit import DNASequence
from dnakit.patterns import scan_restriction_sites

result = scan_restriction_sites(DNASequence("AGAATTCCCGGG"), ["EcoRI"])
hit = result.hits[0]
print(hit.enzyme, hit.top_cut, hit.bottom_cut)
```

- **Example results:**

```text
EcoRI 2 6
```

## 6) PAT-007 CRISPR PAM

- **Meaning:** PAMs are short tags immediately adjacent to CRISPR target sequences, and Cas proteins typically need to recognize PAMs before they can bind and cut nearby DNA.
- **Calculation rules:** Deterministic rules ([FAQ detailed explanation](../../faq.md#pattern-matching-strategy)).
- **Function:** Search for CRISPR guide candidates in the forward and reverse strands according to the specified PAM rules, and return the spacer, PAM, coordinates and direction, which are used for candidate enumeration and do not directly predict editing efficiency or off-target risk.
- **API:** `dnakit.patterns.scan_pam_candidates(value[required], rule[required], guide_length[optional], strand[optional], min_gc[optional], max_gc[optional], exclude_motifs[optional], allow_ambiguous_guides[optional], max_matches[optional], max_scan_length[optional], max_pam_length[optional], max_scan_cells[optional], max_exclude_motifs[optional], max_filter_cells[optional])`, `dnakit.patterns.PAMRule(name[required], pam[required], pam_side[required], guide_length[required], source[optional])`
- **Input:** Required sequence and nuclease name/PAMRule; optional guide length, strand, GC range and exclusion motif.
- **Sample code:**

```python
from dnakit import DNASequence
from dnakit.patterns import scan_pam_candidates

result = scan_pam_candidates(DNASequence("A" * 20 + "TGG"), "SpCas9", strand="forward")
print(result.hits[0].guide_sequence, result.hits[0].pam_sequence)
```

- **Example results:**

```text
AAAAAAAAAAAAAAAAAAAA TGG
```

## 7) PAT-009 palindrome sequence

- **Meaning:** A DNA palindrome is a segment that is identical to its own reverse complement, for example `GAATTC`.
- **Calculation rules:** Deterministic rules ([FAQ detailed explanation](../../faq.md#pattern-matching-strategy)).
- **Function:** Find a palindromic region that is exactly the same as its reverse complement in a sequence, and return the position and length, which can be used to identify possible restriction enzyme sites or symmetrical structural patterns.
- **API:** `dnakit.patterns.find_reverse_complement_palindromes(value[required], min_length[optional], max_length[optional], maximal_per_start[optional], max_comparisons[optional], max_comparison_cells[optional], max_matches[optional])`
- **Input:** Required sequence; optional min/max length, only max hits per starting point, and resource cap.
- **Sample code:**

```python
from dnakit import DNASequence
from dnakit.patterns import find_reverse_complement_palindromes

result = find_reverse_complement_palindromes(
    DNASequence("GAATTC"), min_length=4, max_length=6
)
print(any(hit.sequence == "GAATTC" for hit in result.hits))
```

- **Example results:**

```text
True
```

## 8) PAT-010 Repeat upside down

- **Meaning:** Inverted repeats are two DNA segments that are reverse complementary to each other and may be separated by a sequence. They may fold to form a hairpin structure.
- **Calculation rules:** Deterministic rules ([FAQ detailed explanation](../../faq.md#pattern-matching-strategy)).
- **Function:** Search for an inverted repeat consisting of two reverse complementary sequences and an intermediate loop, and return the coordinates of the two arms and the loop, which is used to screen potential hairpin-related sequence patterns.
- **API:** `dnakit.patterns.find_inverted_repeats(value[required], min_arm_length[optional], max_arm_length[optional], min_loop_length[optional], max_loop_length[optional], max_comparisons[optional], max_comparison_cells[optional], max_matches[optional])`
- **Input:** Required sequence; optional arm/loop length range and comparison, hit upper limit.
- **Sample code:**

```python
from dnakit import DNASequence
from dnakit.patterns import find_inverted_repeats

result = find_inverted_repeats(
    DNASequence("ACGTAAACGT"),
    min_arm_length=4,
    max_arm_length=4,
    min_loop_length=2,
    max_loop_length=2,
)
print(result.hits[0].left_arm, result.hits[0].loop_length)
```

- **Example results:**

```text
ACGT 2
```

## 9) PAT-011 Tandem Repeat

- **Meaning:** Tandem repetition is the same short sequence unit appearing end to end and appearing multiple times in a row. For example, `ATATAT` is `AT` repeated 3 times in a row.
- **Calculation rules:** Deterministic rules ([FAQ detailed explanation](../../faq.md#pattern-matching-strategy)).
- **Function:** Find continuously arranged short repeating units and return the repeating unit, number of times, interval and coverage length, which is used to identify tandem repeats and quantify their size.
- **API:** `dnakit.patterns.find_tandem_repeats(value[required], min_unit_length[optional], max_unit_length[optional], min_repeats[optional], min_repeats_by_unit[optional], overlapping[optional], max_comparisons[optional], max_comparison_cells[optional], max_matches[optional])`
- **Input:** Required sequence; optional unit length, minimum number of repetitions, overlap policy, and resource cap.
- **Sample code:**

```python
from dnakit import DNASequence
from dnakit.patterns import find_tandem_repeats

result = find_tandem_repeats(
    DNASequence("ATATAT"), min_unit_length=1, max_unit_length=3, min_repeats=2
)
print(result.hits[0].unit, result.hits[0].repeat_count)
```

- **Example results:**

```text
AT 3
```

## 10) PAT-012 Microsatellite

- **Meaning:** Microsatellites, also known as STRs, are short tandem repeats with repeating units 1–6 bp in length, such as `CACACACA`.
- **Calculation rules:** Deterministic rules ([FAQ detailed explanation](../../faq.md#pattern-matching-strategy)).
- **Function:** Specifically searches for microsatellites with repeat unit lengths of 1–6 bp, and returns the motif, number of repeats, and coordinates for marking STR candidate regions.
- **API:** `dnakit.patterns.find_microsatellites(value[required], min_repeats_by_unit[optional], overlapping[optional], max_comparisons[optional], max_comparison_cells[optional], max_matches[optional])`
- **Input:** Required sequence; optional minimum number of repetitions per unit length, overlap policy, and resource cap.
- **Sample code:**

```python
from dnakit import DNASequence
from dnakit.patterns import find_microsatellites

result = find_microsatellites(DNASequence("AAAAAACACACA"))
print([(hit.unit, hit.repeat_count) for hit in result.hits])
```

- **Example results:**

```text
[('A', 6), ('CA', 3)]
```
