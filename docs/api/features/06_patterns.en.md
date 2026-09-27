# Sequence function search

Search DNA sequences for patterns such as motifs/PWMs, functional sites, palindromes, inverted repeats, tandem repeats, and microsatellites.

## 1) motif search {#1-pat-001-sim-003-motif-search}

- **Purpose:** Scan the specified motif according to exact, IUPAC or restricted regular rules, and return the coordinates, matching text and chain direction of each hit, which can be used to locate known sequence patterns.
- **Method:** Scan windows using exact, IUPAC or restricted regex matching; PWM scans score each window against the matrix.
- **API:** `dnakit.patterns.scan_motif(value[required], motif[required], mode[optional], name[optional], strand[optional], overlapping[optional], merge_strands[optional], max_matches[optional], max_scan_length[optional], max_pattern_length[optional], max_scan_cells[optional])`, `dnakit.patterns.scan_pwm(value[required], pwm[required], threshold[required], background[optional], pseudocount[optional], strand[optional], max_matches[optional], max_scan_length[optional], max_pwm_length[optional], max_score_cells[optional])`, `dnakit.patterns.PWM(name[required], matrix[required])`
- **Input:** Required sequence and motif/regex/PWM; optional pattern, strand, overlap, background, threshold and hit cap.
- **Output:** Returns `PatternResult` with matched sites and coordinates in `hits`.
- **Example code:**

```python
from dnakit import DNASequence
from dnakit.patterns import scan_motif

result = scan_motif(DNASequence("AAAA"), "AA", strand="forward")
print([(hit.symbol_location.start, hit.symbol_location.end) for hit in result.hits])
```

- **Example output:**

```text
[(0, 2), (1, 3), (2, 4)]
```

## 2) start and stop codons {#2-pat-003-start-and-stop-codons}

- **Purpose:** Search the start and stop codons in the six reading frames of both forward and reverse strands, return the type, reading frame and position, and provide candidate sites for ORF and coding region analysis.
- **Method:** Check triplets in the selected strands and reading frames against the genetic code.
- **API:** `dnakit.patterns.scan_codon_sites(value[required], genetic_code[optional], start_codons[optional], stop_codons[optional], strand[optional], max_matches[optional], max_codon_checks[optional])`
- **Input:** Required sequence; optional genetic code table, custom start/stop set, and strand.
- **Output:** Returns `PatternResult` with matched sites and coordinates in `hits`.
- **Example code:**

```python
from dnakit import DNASequence
from dnakit.patterns import scan_codon_sites

result = scan_codon_sites(DNASequence("ATGAAATAA"), strand="forward")
print([(hit.kind, hit.codon, hit.frame) for hit in result.hits])
```

- **Example output:**

```text
[('start', 'ATG', 1), ('stop', 'TAA', 1), ('stop', 'TGA', 2)]
```

## 3) promoter motif {#3-pat-004-promoter-motif}

- **Purpose:** Scan the built-in consensus sequence or user-provided promoter motif and return candidate positions and chain directions for regular promoter region screening, but does not predict promoter activity.
- **Method:** Match built-in or custom consensus motifs; rule matches do not measure promoter activity.
- **API:** `dnakit.patterns.scan_promoter_motifs(value[required], motifs[optional], strand[optional], max_matches[optional], max_scan_length[optional], max_scan_cells[optional], max_motifs[optional])`
- **Input:** Required sequence; optional motif mapping, strand and resource cap.
- **Output:** Returns `PatternResult` with matched sites and coordinates in `hits`.
- **Example code:**

```python
from dnakit import DNASequence
from dnakit.patterns import scan_promoter_motifs

result = scan_promoter_motifs(DNASequence("GGTATAATCC"), strand="forward")
print([hit.motif_name for hit in result.hits])
```

- **Example output:**

```text
['bacterial_minus_10_consensus']
```

## 4) TF motif {#4-pat-005-tf-motif}

- **Purpose:** Use the position weight matrix provided by the user for window-by-window scoring, and return candidate transcription factor binding sites and scores that exceed the threshold, which are used for motif candidate screening rather than actual binding strength prediction.
- **Method:** Score windows using the PWM, background and pseudocount, retaining scores at or above the threshold.
- **API:** `dnakit.patterns.PWM(name[required], matrix[required])`, `dnakit.patterns.scan_tf_pwm(value[required], tf_name[required], pwm[required], threshold[required], background[optional], pseudocount[optional], strand[optional], max_matches[optional], max_scan_length[optional], max_pwm_length[optional], max_score_cells[optional])`
- **Input:** Required sequence, TF name, PWM and threshold; optional background, pseudocount and strand.
- **Output:** `PWM` stores the weight matrix; `scan_tf_pwm` returns matching sites and scores.
- **Example code:**

```python
from dnakit import DNASequence
from dnakit.patterns import PWM, scan_tf_pwm

pwm = PWM("input", {"A": [4], "C": [0], "G": [0], "T": [0]})
result = scan_tf_pwm(DNASequence("AA"), "TF-X", pwm, threshold=1.0, strand="forward")
print([hit.motif_name for hit in result.hits])
```

- **Example output:**

```text
['TF-X', 'TF-X']
```

## 5) restriction enzyme site {#5-pat-006-restriction-enzyme-site}

- **Purpose:** Search for recognition sequences based on built-in or custom restriction enzyme definitions, and return forward and reverse strand hits and cutting coordinates, which are used to plan enzyme digestion experiments and determine fragment boundaries.
- **Method:** Match recognition sequences and convert enzyme-specific offsets into cut positions on both strands.
- **API:** `dnakit.patterns.scan_restriction_sites(value[required], enzymes[required], max_matches[optional], max_scan_length[optional], max_scan_cells[optional], max_enzymes[optional])`, `dnakit.patterns.RestrictionEnzyme(name[required], recognition_sequence[required], top_cut_offset[required], bottom_cut_offset[required], source[optional])`
- **Input:** Required sequence and enzyme name/definition list; optional scan and hit cap.
- **Output:** Returns `PatternResult` with matched sites and coordinates in `hits`.
- **Example code:**

```python
from dnakit import DNASequence
from dnakit.patterns import scan_restriction_sites

result = scan_restriction_sites(DNASequence("AGAATTCCCGGG"), ["EcoRI"])
hit = result.hits[0]
print(hit.enzyme, hit.top_cut, hit.bottom_cut)
```

- **Example output:**

```text
EcoRI 2 6
```

## 6) CRISPR PAM {#6-pat-007-crispr-pam}

- **Purpose:** Search for CRISPR guide candidates in the forward and reverse strands according to the specified PAM rules, and return the spacer, PAM, coordinates and direction, which are used for candidate enumeration and do not directly predict editing efficiency or off-target risk.
- **Method:** Match PAMs, extract adjacent spacers by orientation, then apply length, GC and excluded-motif filters.
- **API:** `dnakit.patterns.scan_pam_candidates(value[required], rule[required], guide_length[optional], strand[optional], min_gc[optional], max_gc[optional], exclude_motifs[optional], allow_ambiguous_guides[optional], max_matches[optional], max_scan_length[optional], max_pam_length[optional], max_scan_cells[optional], max_exclude_motifs[optional], max_filter_cells[optional])`, `dnakit.patterns.PAMRule(name[required], pam[required], pam_side[required], guide_length[required], source[optional])`
- **Input:** Required sequence and nuclease name/PAMRule; optional guide length, strand, GC range and exclusion motif.
- **Output:** Returns `PatternResult` with matched sites and coordinates in `hits`.
- **Example code:**

```python
from dnakit import DNASequence
from dnakit.patterns import scan_pam_candidates

result = scan_pam_candidates(DNASequence("A" * 20 + "TGG"), "SpCas9", strand="forward")
print(result.hits[0].guide_sequence, result.hits[0].pam_sequence)
```

- **Example output:**

```text
AAAAAAAAAAAAAAAAAAAA TGG
```

## 7) palindrome sequence {#7-pat-009-palindrome-sequence}

- **Purpose:** Find a palindromic region that is exactly the same as its reverse complement in a sequence, and return the position and length, which can be used to identify possible restriction enzyme sites or symmetrical structural patterns.
- **Method:** Enumerate candidate intervals and test equality with their reverse complements.
- **API:** `dnakit.patterns.find_reverse_complement_palindromes(value[required], min_length[optional], max_length[optional], maximal_per_start[optional], max_comparisons[optional], max_comparison_cells[optional], max_matches[optional])`
- **Input:** Required sequence; optional min/max length, only max hits per starting point, and resource cap.
- **Output:** Returns `PatternResult` with matched sites and coordinates in `hits`.
- **Example code:**

```python
from dnakit import DNASequence
from dnakit.patterns import find_reverse_complement_palindromes

result = find_reverse_complement_palindromes(
    DNASequence("GAATTC"), min_length=4, max_length=6
)
print(any(hit.sequence == "GAATTC" for hit in result.hits))
```

- **Example output:**

```text
True
```

## 8) Repeat upside down {#8-pat-010-repeat-upside-down}

- **Purpose:** Search for an inverted repeat consisting of two reverse complementary sequences and an intermediate loop, and return the coordinates of the two arms and the loop, which is used to screen potential hairpin-related sequence patterns.
- **Method:** Enumerate arm and loop lengths and compare both arms for reverse-complement equality.
- **API:** `dnakit.patterns.find_inverted_repeats(value[required], min_arm_length[optional], max_arm_length[optional], min_loop_length[optional], max_loop_length[optional], max_comparisons[optional], max_comparison_cells[optional], max_matches[optional])`
- **Input:** Required sequence; optional arm/loop length range and comparison, hit upper limit.
- **Output:** Returns `PatternResult` with matched sites and coordinates in `hits`.
- **Example code:**

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

- **Example output:**

```text
ACGT 2
```

## 9) Tandem Repeat {#9-pat-011-tandem-repeat}

- **Purpose:** Find continuously arranged short repeating units and return the repeating unit, number of times, interval and coverage length, which is used to identify tandem repeats and quantify their size.
- **Method:** Check consecutive repeats by unit length and report intervals meeting the minimum repeat count.
- **API:** `dnakit.patterns.find_tandem_repeats(value[required], min_unit_length[optional], max_unit_length[optional], min_repeats[optional], min_repeats_by_unit[optional], overlapping[optional], max_comparisons[optional], max_comparison_cells[optional], max_matches[optional])`
- **Input:** Required sequence; optional unit length, minimum number of repetitions, overlap policy, and resource cap.
- **Output:** Returns `PatternResult` with matched sites and coordinates in `hits`.
- **Example code:**

```python
from dnakit import DNASequence
from dnakit.patterns import find_tandem_repeats

result = find_tandem_repeats(
    DNASequence("ATATAT"), min_unit_length=1, max_unit_length=3, min_repeats=2
)
print(result.hits[0].unit, result.hits[0].repeat_count)
```

- **Example output:**

```text
AT 3
```

## 10) Microsatellite {#10-pat-012-microsatellite}

- **Purpose:** Specifically searches for microsatellites with repeat unit lengths of 1–6 bp, and returns the motif, number of repeats, and coordinates for marking STR candidate regions.
- **Method:** Scan tandem repeats with unit lengths restricted to 1–6 bp.
- **API:** `dnakit.patterns.find_microsatellites(value[required], min_repeats_by_unit[optional], overlapping[optional], max_comparisons[optional], max_comparison_cells[optional], max_matches[optional])`
- **Input:** Required sequence; optional minimum number of repetitions per unit length, overlap policy, and resource cap.
- **Output:** Returns `PatternResult` with matched sites and coordinates in `hits`.
- **Example code:**

```python
from dnakit import DNASequence
from dnakit.patterns import find_microsatellites

result = find_microsatellites(DNASequence("AAAAAACACACA"))
print([(hit.unit, hit.repeat_count) for hit in result.hits])
```

- **Example output:**

```text
[('A', 6), ('CA', 3)]
```

[FAQ: details](../../faq/tutorial_basics.md#pattern-terms)

## 11) Build a motif and scan DNA {#bio-49}

- **Purpose:** Build a motif model from binding sites and locate matches above a threshold.
- **Method:** Compute counts/PWM/PSSM, derive a conservative FPR threshold from the background model, and scan either strand or circular sequences.
- **API:** `dnakit.bioinformatics.motifs.build_motif(sites=None, counts=None, pseudocount=0.5)`; `dnakit.bioinformatics.motifs.fpr_threshold(model, fpr=0.0001)`; `dnakit.bioinformatics.motifs.scan_motif_model(sequence, model, threshold=..., both_strands=True, circular=False)`
- **Input:** Equal-width binding sites or a counts matrix, plus a target sequence; requires `bio-analysis`.
- **Output:** `build_motif` returns a `MotifModel`; scanning returns hits with positions and scores.
- **Example code:**

```python
from dnakit.bioinformatics.motifs import build_motif, fpr_threshold, scan_motif_model

model = build_motif(["ACG", "ACG", "ATG"], pseudocount=0.1)
threshold = fpr_threshold(model, fpr=0.1, precision=100)
hits = scan_motif_model("TTACGTT", model, threshold=threshold, both_strands=False)
print(model.consensus)
print([(hit["start"], hit["end"]) for hit in hits])
```

- **Example output:**

```text
ACG
[(2, 5)]
```

Dependencies and details: [FAQ](../../faq/bioinformatics.md#optimization).
