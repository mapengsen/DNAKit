# 72 Variants: offline varcode, hgvs and scikit-allel functions

This page covers the varcode, hgvs and scikit-allel functions wrapped in batch G6a. Every function reads local inputs only; none needs a reference genome, a pyensembl cache, seqrepo, UTA, PostgreSQL or network access.

The wrappers live in `dnakit.variants`, `dnakit.variants` and `dnakit.variants`; the examples import from these modules.

Dependencies are optional packages: `bio-varcode`, `bio-hgvs`, `bio-scikit-allel` and `bio-analysis` (numpy). A missing package raises `BackendUnavailableError` with an install hint.

## Overview

| Package | Upstream public entries | Wrapped in this batch | Already wrapped by dnakit | Family variants | Duplicates of dnakit functions |
| --- | ---: | ---: | ---: | ---: | ---: |
| varcode | 333 | 18 | 1 | 2 | 6 |
| hgvs | 31 | 6 | 0 | 2 | 1 |
| scikit-allel | 634 | 20 | 11 | 24 | 6 |

Other entries record their reason in the ledger: they need pyensembl or Ensembl data, seqrepo or UTA, duplicate an existing dnakit function, are deprecated, are command-line entry points, or fall beyond the wrapper cap.

## 1) Enumerate cryptic exon candidates

- **Purpose:** Find acceptor-donor pairs in a local sequence that could form a new internal exon.
- **Method:** Score every 9 bp donor and 4 bp acceptor window with varcode's consensus model, keep sites at or above the thresholds, pair them within the intron length range, and sort by summed score, start and end.
- **API:** `dnakit.variants.enumerate_cryptic_exon_candidates(contig: str, sequence: str, sequence_start: int = 1, *, donor_threshold: float = 0.7, acceptor_threshold: float = 0.7, min_intron_length: int = 30, max_intron_length: int = 10000)`
- **Input:** Contig, sequence (A/C/G/T/N only), 1-based sequence start; optional thresholds and intron length range.
- **Output:** `candidates` list with contig, interval_start, interval_end, donor_score and acceptor_score.
- **Example code:**

```python
from dnakit.variants import enumerate_cryptic_exon_candidates

sequence = "TTTT" + "CAGG" + "C" * 30 + "AAG" + "GTAAGT" + "TTTT"
result = enumerate_cryptic_exon_candidates("chr1", sequence, 1000, min_intron_length=10)
print(result["count"], result["candidates"][0])
```

- **Example output:**

```text
1 {'contig': 'chr1', 'interval_start': 1007, 'interval_end': 1041, 'donor_score': 1.0, 'acceptor_score': 1.0}
```

## 2) Score a splice-site window

- **Purpose:** Give a consensus match score to a 9 bp donor or 4 bp acceptor window.
- **Method:** Check the GT (donor) or AG (acceptor) anchor bases; a failed anchor scores 0, otherwise the score is the fraction of bases matching the consensus.
- **API:** `dnakit.variants.score_splice_site(window: str, *, site: str = 'donor')`
- **Input:** Window sequence; `site` is `donor` (9 bp) or `acceptor` (4 bp).
- **Output:** `score` between 0 and 1 and the window length.
- **Example code:**

```python
from dnakit.variants import score_splice_site

print(score_splice_site("CAGGTAAGT", site="donor")["score"])
print(score_splice_site("CAGG", site="acceptor")["score"])
print(score_splice_site("CAGATAAGT", site="donor")["score"])
```

- **Example output:**

```text
1.0
1.0
0.0
```

## 3) Parse a symbolic SV ALT

- **Purpose:** Turn `<DEL>`, `<DUP>`, `<INV>`, `<CNV>` or breakend ALTs into structural-variant fields.
- **Method:** Read END, SVTYPE and MATEID from INFO as in the VCF convention; an END not greater than POS is rejected.
- **API:** `dnakit.variants.parse_symbolic_sv_alt(contig: str, start: int, ref: str, alt: str, *, info: Mapping[str, Any] | None = None)`
- **Input:** Contig, POS, REF, ALT and an optional INFO mapping.
- **Output:** `sv_type`, `start`, `end`, `affected_start` and `affected_end`; breakend records also carry `mate_contig` and `mate_start`.
- **Example code:**

```python
from dnakit.variants import parse_symbolic_sv_alt

record = parse_symbolic_sv_alt("chr1", 100, "N", "<DEL>", info={"END": 500})
print(record["sv_type"], record["start"], record["end"], record["affected_start"], record["affected_end"])
```

- **Example output:**

```text
DEL 100 500 101 500
```

## 4) Read which side a breakend ALT keeps

- **Purpose:** Tell which side of each breakpoint a VCF breakend ALT keeps.
- **Method:** Match the bracket and prefix or suffix form; forms such as `t]p]` map to left or right.
- **API:** `dnakit.variants.breakend_alt_sides(alt: str)`
- **Input:** VCF breakend ALT string.
- **Output:** `valid`, `this_side` and `mate_side`; a non-breakend ALT gives `valid=False`.
- **Example code:**

```python
from dnakit.variants import breakend_alt_sides

sides = breakend_alt_sides("A]chr1:500]")
print(sides["valid"], sides["this_side"], sides["mate_side"])
```

- **Example output:**

```text
True left left
```

## 5) Read the mate position of a breakend ALT

- **Purpose:** Extract the mate contig, position and bracket orientation from a breakend ALT.
- **Method:** Parse the ALT with the VCF breakend pattern; an invalid ALT returns `valid=False`.
- **API:** `dnakit.variants.breakend_mate_position(alt: str)`
- **Input:** VCF breakend ALT string.
- **Output:** `mate_contig`, `mate_position` and `mate_orientation`.
- **Example code:**

```python
from dnakit.variants import breakend_mate_position

mate = breakend_mate_position("A]chr1:500]")
print(mate["mate_contig"], mate["mate_position"], mate["mate_orientation"])
```

- **Example output:**

```text
chr1 500 ]]
```

## 6) Read bases inserted by a breakend ALT

- **Purpose:** Give the bases a breakend ALT inserts beside a one-base REF anchor.
- **Method:** Prefix form keeps bases after the anchor, suffix form bases before it; a multi-base REF or a mismatched anchor is unknown.
- **API:** `dnakit.variants.breakend_insert_bases(ref: str, alt: str)`
- **Input:** One-base REF and a breakend ALT.
- **Output:** `known` and `inserted_bases`; `inserted_bases` is None when unknown.
- **Example code:**

```python
from dnakit.variants import breakend_insert_bases

inserted = breakend_insert_bases("A", "AGT[chr1:500[")
print(inserted["known"], inserted["inserted_bases"])
```

- **Example output:**

```text
True GT
```

## 7) Type a pair of breakend records

- **Purpose:** Classify two paired breakend records as DEL, DUP or INV.
- **Method:** The records must point at each other on one contig with the same SVTYPE; varcode's paired-breakend typing decides.
- **API:** `dnakit.variants.classify_breakend_pair(first: Mapping[str, Any], second: Mapping[str, Any])`
- **Input:** Two breakend records, each with contig, start, ref, alt and optional info.
- **Output:** `typed`, `sv_type`, `start` and `end` (`start` is the base before the event).
- **Example code:**

```python
from dnakit.variants import classify_breakend_pair

first = {"contig": "chr1", "start": 100, "ref": "A", "alt": "A]chr1:500]",
         "info": {"SVTYPE": "DEL", "MATEID": "b2"}}
second = {"contig": "chr1", "start": 500, "ref": "T", "alt": "[chr1:100[T",
          "info": {"SVTYPE": "DEL", "MATEID": "b1"}}
typed = classify_breakend_pair(first, second)
print(typed["sv_type"], typed["start"], typed["end"])
```

- **Example output:**

```text
DEL 100 499
```

## 8) Compare structural-variant call breakpoints

- **Purpose:** Group structural-variant calls whose breakpoints are close, keeping every input record.
- **Method:** Bucket calls by normalized breakend topology, then greedily group within a bucket with complete linkage; calls farther apart than `max_distance` stay separate.
- **API:** `dnakit.variants.compare_sv_calls(rows: Sequence[Mapping[str, Any]], *, max_distance: int = 100)`
- **Input:** List of call records with call_id, caller, sample, build (GRCh37 or GRCh38), chrom, pos, ref and alt.
- **Output:** `result['groups']`, each with call_ids, callers, status and related fields.
- **Example code:**

```python
from dnakit.variants import compare_sv_calls

rows = [
    {"call_id": "a", "caller": "x", "sample": "s", "build": "GRCh38", "chrom": "chr1",
     "pos": 100, "ref": "A", "alt": "<DEL>", "info": {"END": 500}},
    {"call_id": "b", "caller": "y", "sample": "s", "build": "GRCh38", "chrom": "chr1",
     "pos": 100, "ref": "A", "alt": "<DEL>", "info": {"END": 500}},
    {"call_id": "c", "caller": "x", "sample": "s", "build": "GRCh38", "chrom": "chr1",
     "pos": 5000, "ref": "A", "alt": "<DEL>", "info": {"END": 5400}},
]
result = compare_sv_calls(rows, max_distance=100)
print(sorted(group["call_ids"] for group in result["result"]["groups"]))
```

- **Example output:**

```text
[['a', 'b'], ['c']]
```

## 9) Trim shared allele flanks

- **Purpose:** Remove the prefix and suffix shared by REF and ALT to get the minimal allele form.
- **Method:** Strip the common prefix first, then the common suffix of what remains, as varcode does before classifying.
- **API:** `dnakit.variants.trim_shared_allele_flanks(ref: str, alt: str)`
- **Input:** REF and ALT strings (either may be empty).
- **Output:** `trimmed_ref`, `trimmed_alt`, `prefix` and `suffix`.
- **Example code:**

```python
from dnakit.variants import trim_shared_allele_flanks

result = trim_shared_allele_flanks("ACGT", "ACTT")
print(result["trimmed_ref"], result["trimmed_alt"], result["prefix"], result["suffix"])
```

- **Example output:**

```text
G T AC T
```

## 10) Parse a VCF genotype string

- **Purpose:** Parse GT strings such as `0|1` or `./.` into allele indices and phase.
- **Method:** Check the VCF GT pattern first, then parse with varcode; missing alleles are None.
- **API:** `dnakit.variants.parse_vcf_genotype(genotype: str)`
- **Input:** GT string.
- **Output:** `alleles` list, `phased` and `ploidy`.
- **Example code:**

```python
from dnakit.variants import parse_vcf_genotype

genotype = parse_vcf_genotype("1|.")
print(genotype["alleles"], genotype["phased"], genotype["ploidy"])
```

- **Example output:**

```text
[1, None] True 2
```

## 11) Normalize an allele string

- **Purpose:** Uppercase an allele and record gap characters such as `.` or `-` as an empty allele.
- **Method:** Call varcode's nucleotide normalizer; by default only standard bases and gap characters are accepted.
- **API:** `dnakit.variants.normalize_nucleotide_allele(allele: str, *, allow_extended_nucleotides: bool = False, empty_chars: str = '.-', treat_nan_as_empty: bool = True)`
- **Input:** Allele string; optional `allow_extended_nucleotides` and `empty_chars`.
- **Output:** `normalized` and `is_empty`.
- **Example code:**

```python
from dnakit.variants import normalize_nucleotide_allele

print(repr(normalize_nucleotide_allele("acgt")["normalized"]))
print(normalize_nucleotide_allele(".")["is_empty"])
```

- **Example output:**

```text
'ACGT'
True
```

## 12) Test whether a variant overlaps an interval

- **Purpose:** Check whether a variant of n reference bases meets a 1-based inclusive interval.
- **Method:** The variant covers start to start+n-1; an insertion (n=0) overlaps only when the interval contains the insertion point.
- **API:** `dnakit.variants.variant_overlaps_interval(variant_start: int, n_ref_bases: int, interval_start: int, interval_end: int)`
- **Input:** Variant start, reference base count, interval start and end (1-based, inclusive).
- **Output:** `overlaps` and the variant end `variant_end`.
- **Example code:**

```python
from dnakit.variants import variant_overlaps_interval

print(variant_overlaps_interval(10, 3, 12, 20)["overlaps"])
print(variant_overlaps_interval(10, 3, 13, 20)["overlaps"])
```

- **Example output:**

```text
True
False
```

## 13) Four-gamete test for cis or trans

- **Purpose:** Decide whether two variants lie on the same chromosome (cis) or on different ones (trans) from read-fragment counts.
- **Method:** Call varcode's four-gamete test on the AB, Ab, aB and ab counts; `max_p_value` sets significance and sparse evidence returns unknown.
- **API:** `dnakit.variants.four_gamete_phase_test(both: int, first: int, second: int, neither: int, *, min_fragments: int = 2, error_rate: float = 0.01, max_p_value: float = 0.05)`
- **Input:** Four non-negative counts (both, first, second, neither) and test parameters.
- **Output:** `relation` (cis, trans or unknown) and `in_cis`.
- **Example code:**

```python
from dnakit.variants import four_gamete_phase_test

cis = four_gamete_phase_test(10, 1, 2, 10)
print(cis["relation"], cis["in_cis"])
```

- **Example output:**

```text
cis True
```

## 14) Detect germline allele overlap with a somatic variant

- **Purpose:** Find germline variants that share the contig, position, REF and ALT of a somatic variant.
- **Method:** Call varcode's overlap test for each germline record; genotype, allele balance and LOH are not considered.
- **API:** `dnakit.variants.detect_germline_allele_overlap(somatic: Mapping[str, Any], germline: Sequence[Mapping[str, Any]])`
- **Input:** A somatic variant and a list of germline variants, each with contig, start, ref and alt.
- **Output:** `overlap`, `matched_indices` and `germline_count`.
- **Example code:**

```python
from dnakit.variants import detect_germline_allele_overlap

somatic = {"contig": "chr1", "start": 100, "ref": "A", "alt": "G"}
germline = [{"contig": "chr1", "start": 100, "ref": "A", "alt": "G"},
            {"contig": "chr1", "start": 200, "ref": "C", "alt": "T"}]
result = detect_germline_allele_overlap(somatic, germline)
print(result["overlap"], result["matched_indices"])
```

- **Example output:**

```text
True [0]
```

## 15) Partition germline variants by phase

- **Purpose:** Split germline variants into groups in cis or trans with a somatic variant, linked groups, or unphased.
- **Method:** Homozygous variants are always cis; the rest propagate through supplied in-cis evidence with union-find, and variants without evidence stay in `unphased`.
- **API:** `dnakit.variants.partition_germline_by_phase(somatic: Mapping[str, Any], germline: Sequence[Mapping[str, Any]], *, homozygous: Sequence[Mapping[str, Any]] | None = None, cis_variants: Sequence[Mapping[str, Any]] | None = None, in_cis_evidence: Sequence[Mapping[str, Any]] | None = None)`
- **Input:** Somatic variant, germline variants, and optional homozygous variants and `in_cis_evidence` (each with first, second and in_cis).
- **Output:** `cis`, `trans`, `linked` and `unphased`.
- **Example code:**

```python
from dnakit.variants import partition_germline_by_phase

somatic = {"contig": "chr1", "start": 100, "ref": "A", "alt": "G"}
germline = {"contig": "chr1", "start": 200, "ref": "C", "alt": "T"}
evidence = [{"first": somatic, "second": germline, "in_cis": False}]
partition = partition_germline_by_phase(somatic, [germline], in_cis_evidence=evidence)
print(partition["cis"], partition["trans"], partition["unphased"])
```

- **Example output:**

```text
[] [{'contig': 'chr1', 'start': 200, 'ref': 'C', 'alt': 'T'}] []
```

## 16) Enumerate phase hypotheses

- **Purpose:** List the phase configurations a somatic variant can take against germline variants.
- **Method:** Partition the germline variants, then build one hypothesis per configuration given the in-cis evidence; more than `max_hypotheses` raises an error.
- **API:** `dnakit.variants.enumerate_phase_hypotheses(somatic: Mapping[str, Any], germline: Sequence[Mapping[str, Any]], *, in_cis_evidence: Sequence[Mapping[str, Any]] | None = None, max_hypotheses: int = 8)`
- **Input:** Somatic variant, germline variants, optional `in_cis_evidence`, and `max_hypotheses` (1 to 64).
- **Output:** `count` and `hypotheses`, each with phase_state, haplotype, cis and trans.
- **Example code:**

```python
from dnakit.variants import enumerate_phase_hypotheses

somatic = {"contig": "chr1", "start": 100, "ref": "A", "alt": "G"}
germline = {"contig": "chr1", "start": 200, "ref": "C", "alt": "T"}
evidence = [{"first": somatic, "second": germline, "in_cis": True}]
result = enumerate_phase_hypotheses(somatic, [germline], in_cis_evidence=evidence)
print(result["count"], result["hypotheses"][0]["phase_state"])
```

- **Example output:**

```text
1 phased
```

## 17) Load a local MAF file

- **Purpose:** Read a TCGA-style MAF file into variant records with contig, position, REF and ALT.
- **Method:** Call varcode's MAF reader; records are sorted by contig and position and keep their reference_name.
- **API:** `dnakit.variants.load_maf_variants(path: str)`
- **Input:** Path to a local MAF file with the 19 standard columns varcode requires.
- **Output:** `reference`, `count` and `variants`.
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.variants import load_maf_variants

COLUMNS = ["Hugo_Symbol", "Entrez_Gene_Id", "Center", "NCBI_Build", "Chromosome", "Start_Position",
           "End_Position", "Strand", "Variant_Classification", "Variant_Type", "Reference_Allele",
           "Tumor_Seq_Allele1", "Tumor_Seq_Allele2", "dbSNP_RS", "dbSNP_Val_Status",
           "Tumor_Sample_Barcode", "Matched_Norm_Sample_Barcode", "Match_Norm_Seq_Allele1",
           "Match_Norm_Seq_Allele2"]
ROWS = [
    {"Hugo_Symbol": "TP53", "NCBI_Build": "GRCh38", "Chromosome": "17", "Start_Position": "7676154",
     "End_Position": "7676154", "Strand": "+", "Variant_Type": "SNP", "Reference_Allele": "C",
     "Tumor_Seq_Allele1": "C", "Tumor_Seq_Allele2": "T", "Tumor_Sample_Barcode": "S1"},
    {"Hugo_Symbol": "KRAS", "NCBI_Build": "GRCh38", "Chromosome": "12", "Start_Position": "25245350",
     "End_Position": "25245350", "Strand": "+", "Variant_Type": "SNP", "Reference_Allele": "C",
     "Tumor_Seq_Allele1": "C", "Tumor_Seq_Allele2": "A", "Tumor_Sample_Barcode": "S1"},
]
with tempfile.TemporaryDirectory() as folder:
    path = Path(folder) / "sample.maf"
    lines = ["\t".join(COLUMNS)]
    lines += ["\t".join(row.get(column, "") for column in COLUMNS) for row in ROWS]
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    result = load_maf_variants(str(path))
print(result["reference"], result["count"], result["variants"][0]["contig"], result["variants"][0]["start"])
```

- **Example output:**

```text
GRCh38 2 12 25245350
```

## 18) Check sample identity in a VCF

- **Purpose:** Load sample fingerprints from a local VCF and report germline concordance between sample pairs.
- **Method:** Call varcode's fingerprint loader and identity check; the result is a screening heuristic, not an identity probability.
- **API:** `dnakit.variants.check_vcf_sample_identity(path: str, *, kind: str = 'germline', assembly: str | None = None, config: Mapping[str, Any] | None = None)`
- **Input:** Local VCF path; optional `kind` (germline, somatic or mixed), `assembly`, and `config` overrides.
- **Output:** `sample_count` and `report`, which holds germline statistics for each pair.
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.variants import check_vcf_sample_identity

VCF = (
    "##fileformat=VCFv4.2\n"
    "##reference=GRCh38\n"
    "##contig=<ID=chr1,length=248956422>\n"
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tA\tB\n"
    "chr1\t100\t.\tA\tG\t50\tPASS\t.\tGT:DP\t0/1:30\t0/1:30\n"
    "chr1\t200\t.\tC\tT\t50\tPASS\t.\tGT:DP\t1/1:30\t1/1:30\n"
)
with tempfile.TemporaryDirectory() as folder:
    path = Path(folder) / "pair.vcf"
    path.write_text(VCF, encoding="utf-8")
    result = check_vcf_sample_identity(str(path))
germline = result["report"]["pairs"][0]["germline"]
print(result["sample_count"], germline["shared_snps"], germline["exact_matches"])
```

- **Example output:**

```text
2 2 2
```

## 19) Windowed Tajima's D

- **Purpose:** Compute Tajima's D in base-pair windows along one contig to look for selection signals.
- **Method:** Call scikit-allel's windowed_tajima_d; windows are placed by position and need at least `min_sites` sites.
- **API:** `dnakit.variants.calculate_windowed_tajimas_d(positions: Sequence[int], allele_counts: Sequence[Sequence[int]], *, window_size: int, step: int | None = None, min_sites: int = 3)`
- **Input:** 1-based positions in ascending order and allele counts per site (variants x alleles).
- **Output:** `windows`, each with start, stop, sites and tajima_d.
- **Example code:**

```python
from dnakit.variants import calculate_windowed_tajimas_d

positions = [100, 200, 300, 400, 500, 600, 700, 800, 900, 1000]
counts = [[3, 1], [2, 2], [4, 0], [1, 3], [0, 4], [2, 2], [3, 1], [1, 3], [2, 2], [4, 0]]
result = calculate_windowed_tajimas_d(positions, counts, window_size=500)
for window in result["windows"]:
    print(window["start"], window["stop"], window["sites"], round(window["tajima_d"], 6))
```

- **Example output:**

```text
100 599 5 0.167656
600 1000 5 0.650102
```

## 20) Windowed nucleotide diversity

- **Purpose:** Compute nucleotide diversity pi per window, reporting the callable bases used.
- **Method:** Call scikit-allel's windowed_diversity; pi is the summed mean pairwise differences divided by callable bases.
- **API:** `dnakit.variants.calculate_windowed_nucleotide_diversity(positions: Sequence[int], allele_counts: Sequence[Sequence[int]], *, window_size: int, step: int | None = None)`
- **Input:** Positions, allele counts and window size in base pairs.
- **Output:** `windows`, each with start, stop, callable_bases, sites and pi.
- **Example code:**

```python
from dnakit.variants import calculate_windowed_nucleotide_diversity

positions = [100, 200, 300, 400, 500, 600, 700, 800, 900, 1000]
counts = [[3, 1], [2, 2], [4, 0], [1, 3], [0, 4], [2, 2], [3, 1], [1, 3], [2, 2], [4, 0]]
result = calculate_windowed_nucleotide_diversity(positions, counts, window_size=500)
for window in result["windows"]:
    print(window["start"], window["stop"], window["callable_bases"], round(window["pi"], 8))
```

- **Example output:**

```text
100 599 500 0.00333333
600 1000 401 0.00581879
```

## 21) Mean pairwise difference per site

- **Purpose:** Compute the mean number of pairwise allele differences at each site.
- **Method:** From allele counts as 2·n_a·n_b / (n·(n−1)) for a biallelic site; scikit-allel computes the values.
- **API:** `dnakit.variants.calculate_mean_pairwise_difference(allele_counts: Sequence[Sequence[int]])`
- **Input:** Allele counts per site (variants x alleles).
- **Output:** `values`, one per input site.
- **Example code:**

```python
from dnakit.variants import calculate_mean_pairwise_difference

result = calculate_mean_pairwise_difference([[3, 1], [2, 2], [4, 0]])
print([round(value, 6) for value in result["values"]])
```

- **Example output:**

```text
[0.5, 0.666667, 0.0]
```

## 22) Hudson Fst per site

- **Purpose:** Give Hudson's Fst numerator and denominator per site and the overall Fst.
- **Method:** Call scikit-allel's hudson_fst; the overall value is the sum of numerators over the sum of denominators.
- **API:** `dnakit.variants.calculate_hudson_fst_per_site(first_counts: Sequence[Sequence[int]], second_counts: Sequence[Sequence[int]])`
- **Input:** Allele counts for two populations with matching site order.
- **Output:** `numerator`, `denominator` and `fst`.
- **Example code:**

```python
from dnakit.variants import calculate_hudson_fst_per_site

first = [[4, 0], [3, 1], [2, 2], [0, 4]]
second = [[0, 4], [1, 3], [2, 2], [4, 0]]
result = calculate_hudson_fst_per_site(first, second)
print([round(value, 6) for value in result["numerator"]], round(result["fst"], 6))
```

- **Example output:**

```text
[1.0, 0.125, -0.166667, 1.0] 0.626667
```

## 23) Patterson's D statistic

- **Purpose:** Compute the ABBA-BABA D statistic from allele counts of four populations.
- **Method:** Call scikit-allel's patterson_d; the numerator is (a−b)(c−d) and the denominator the variance term; D is the summed numerator over the summed denominator.
- **API:** `dnakit.variants.calculate_patterson_d(first_counts: Sequence[Sequence[int]], second_counts: Sequence[Sequence[int]], third_counts: Sequence[Sequence[int]], fourth_counts: Sequence[Sequence[int]])`
- **Input:** Allele counts for populations A, B, C and D with matching site order.
- **Output:** `numerator`, `denominator` and `d`.
- **Example code:**

```python
from dnakit.variants import calculate_patterson_d

result = calculate_patterson_d([[4, 0]], [[0, 4]], [[4, 0]], [[0, 4]])
print(result["numerator"], result["denominator"], result["d"])
```

- **Example output:**

```text
[1.0] [1.0] 1.0
```

## 24) Expected heterozygosity

- **Purpose:** Compute expected heterozygosity 1 − Σp² from allele frequencies at each site.
- **Method:** Diploid Hardy-Weinberg expectation; sites without valid frequencies give None.
- **API:** `dnakit.variants.calculate_expected_heterozygosity(allele_frequencies: Sequence[Sequence[float]])`
- **Input:** Allele frequency rows per site; each row sums to 1.
- **Output:** `values` and `ploidy`, which is fixed at 2.
- **Example code:**

```python
from dnakit.variants import calculate_expected_heterozygosity

result = calculate_expected_heterozygosity([[0.5, 0.5], [0.75, 0.25], [1.0, 0.0]])
print(result["values"])
```

- **Example output:**

```text
[0.5, 0.375, 0.0]
```

## 25) Inbreeding coefficient

- **Purpose:** Compute the inbreeding coefficient F = 1 − Ho / He for each site.
- **Method:** Observed heterozygosity Ho and expected He from diploid genotypes; sites with He equal to zero give None.
- **API:** `dnakit.variants.calculate_inbreeding_coefficient(genotypes: Sequence[Sequence[Sequence[int]]])`
- **Input:** Diploid genotype array (variants x samples x 2); missing calls are −1.
- **Output:** `values` per site.
- **Example code:**

```python
from dnakit.variants import calculate_inbreeding_coefficient

genotypes = [[[0, 0], [0, 1]], [[1, 1], [0, 1]], [[0, 1], [1, 1]], [[0, 0], [0, 0]]]
print([None if value is None else round(value, 6) for value in calculate_inbreeding_coefficient(genotypes)["values"]])
```

- **Example output:**

```text
[-0.333333, -0.333333, -0.333333, None]
```

## 26) Rogers-Huff r

- **Purpose:** Compute Rogers-Huff r, the correlation of genotype dosages, for every pair of sites.
- **Method:** Pearson correlation of the dosage vectors; a site with zero variance gives None.
- **API:** `dnakit.variants.calculate_rogers_huff_r(alt_counts: Sequence[Sequence[int]])`
- **Input:** Allele dosage matrix (variants x samples) with values 0, 1, 2 or −1 for missing; at most 200 sites.
- **Output:** `pairs`, each with variant_a, variant_b and r.
- **Example code:**

```python
from dnakit.variants import calculate_rogers_huff_r

result = calculate_rogers_huff_r([[0, 1], [2, 1], [1, 2], [0, 0]])
for pair in result["pairs"]:
    print(pair["variant_a"], pair["variant_b"], pair["r"])
```

- **Example output:**

```text
0 1 -1.0
0 2 1.0
0 3 None
1 2 -1.0
1 3 None
2 3 None
```

## 27) Select sites in approximate linkage equilibrium

- **Purpose:** Prune sites with a sliding-window r-squared filter, keeping sites close to linkage equilibrium.
- **Method:** Call scikit-allel's locate_unlinked; the window and step count sites, and sites above the threshold are removed.
- **API:** `dnakit.variants.select_unlinked_variants(alt_counts: Sequence[Sequence[int]], *, size: int = 100, step: int = 20, threshold: float = 0.1)`
- **Input:** Dosage matrix (variants x samples) plus window size, step and threshold.
- **Output:** `kept_indices` in ascending order and `kept_count`.
- **Example code:**

```python
import numpy as np

from dnakit.variants import select_unlinked_variants

rng = np.random.default_rng(1)
dosage = rng.integers(0, 3, size=(30, 8))
result = select_unlinked_variants(dosage.tolist(), size=10, step=5, threshold=0.1)
print(result["kept_indices"])
```

- **Example output:**

```text
[0, 1, 5, 12, 14, 15, 20, 27, 28]
```

## 28) Pairwise sample distances

- **Purpose:** Compute a distance matrix between samples from site dosages.
- **Method:** Call scikit-allel's pairwise_distance for the condensed form and expand it to a symmetric square; missing calls are rejected.
- **API:** `dnakit.variants.calculate_pairwise_distance(alt_counts: Sequence[Sequence[int]], *, metric: str = 'euclidean')`
- **Input:** Non-negative dosage matrix (variants x samples); optional metric.
- **Output:** `matrix` (samples x samples) and `metric`.
- **Example code:**

```python
from dnakit.variants import calculate_pairwise_distance

result = calculate_pairwise_distance([[0, 1, 2], [1, 1, 2], [0, 0, 1]])
print([[round(value, 6) for value in row] for row in result["matrix"]])
```

- **Example output:**

```text
[[0.0, 1.0, 2.44949], [1.0, 0.0, 1.732051], [2.44949, 1.732051, 0.0]]
```

## 29) Unstandardized iHS

- **Purpose:** Compute the unstandardized integrated haplotype score (iHS) for each site in one population.
- **Method:** Call scikit-allel's ihs, integrating EHH decay; sites with allele frequency below `min_maf` are not scored.
- **API:** `dnakit.variants.calculate_ihs(haplotypes: Sequence[Sequence[int]], positions: Sequence[int], *, min_ehh: float = 0.05, min_maf: float = 0.05, gap_scale: int = 20000, max_gap: int = 200000, include_edges: bool = False)`
- **Input:** Haplotype matrix (variants x haplotypes, values 0 or 1), positions and EHH parameters.
- **Output:** `ihs` list, one value per site.
- **Example code:**

```python
import numpy as np

from dnakit.variants import calculate_ihs

rng = np.random.default_rng(3)
haplotypes = rng.integers(0, 2, size=(40, 12))
positions = list(range(1000, 1000 + 100 * 40, 100))
values = calculate_ihs(haplotypes.tolist(), positions)["ihs"]
print(len(values), [round(value, 6) for value in values[5:8]])
```

- **Example output:**

```text
40 [0.143101, -0.169076, 0.077962]
```

## 30) Unstandardized XP-EHH

- **Purpose:** Compare haplotypes from two populations and give the cross-population EHH score per site.
- **Method:** Call scikit-allel's xpehh; both haplotype sets must describe the same sites.
- **API:** `dnakit.variants.calculate_xpehh(haplotypes_a: Sequence[Sequence[int]], haplotypes_b: Sequence[Sequence[int]], positions: Sequence[int], *, min_ehh: float = 0.05, gap_scale: int = 20000, max_gap: int = 200000, include_edges: bool = False)`
- **Input:** Two haplotype matrices, positions and EHH parameters.
- **Output:** `xpehh` list.
- **Example code:**

```python
import numpy as np

from dnakit.variants import calculate_xpehh

rng = np.random.default_rng(4)
first = rng.integers(0, 2, size=(40, 12))
second = rng.integers(0, 2, size=(40, 12))
positions = list(range(1000, 1000 + 100 * 40, 100))
values = calculate_xpehh(first.tolist(), second.tolist(), positions)["xpehh"]
print(len(values), [round(value, 6) for value in values[5:8]])
```

- **Example output:**

```text
40 [0.157004, 0.055152, 0.196028]
```

## 31) Unstandardized nSL

- **Purpose:** Compute the number of segregating sites by length (nSL) for each site.
- **Method:** Call scikit-allel's nsl, comparing shared haplotype prefix lengths between haplotypes at each site.
- **API:** `dnakit.variants.calculate_nsl(haplotypes: Sequence[Sequence[int]])`
- **Input:** Haplotype matrix (variants x haplotypes, values 0 or 1).
- **Output:** `nsl` list.
- **Example code:**

```python
import numpy as np

from dnakit.variants import calculate_nsl

rng = np.random.default_rng(6)
haplotypes = rng.integers(0, 2, size=(40, 12))
values = calculate_nsl(haplotypes.tolist())["nsl"]
print(len(values), [round(value, 6) for value in values[:3]])
```

- **Example output:**

```text
40 [-0.02353, 0.092288, 0.059423]
```

## 32) Garud H statistics

- **Purpose:** Compute H1, H12, H123 and H2/H1 from haplotype frequencies to detect soft sweeps.
- **Method:** Call scikit-allel's garud_h: H1 is the sum of squared haplotype frequencies; H12 and H123 merge the two or three most frequent haplotypes.
- **API:** `dnakit.variants.calculate_garud_h(haplotypes: Sequence[Sequence[int]])`
- **Input:** Haplotype matrix (variants x haplotypes, values 0 or 1).
- **Output:** `h1`, `h12`, `h123` and `h2_h1`.
- **Example code:**

```python
import numpy as np

from dnakit.variants import calculate_garud_h

rng = np.random.default_rng(7)
haplotypes = rng.integers(0, 2, size=(40, 12))
result = calculate_garud_h(haplotypes.tolist())
print(round(result["h1"], 6), round(result["h12"], 6), round(result["h2_h1"], 6))
```

- **Example output:**

```text
0.083333 0.097222 0.916667
```

## 33) EHH decay curve

- **Purpose:** Compute how extended haplotype homozygosity (EHH) decays moving away from the first site.
- **Method:** Call scikit-allel's ehh_decay; `truncate=True` drops trailing zeros.
- **API:** `dnakit.variants.calculate_ehh_decay(haplotypes: Sequence[Sequence[int]], *, truncate: bool = False)`
- **Input:** Haplotype matrix (variants x haplotypes, values 0 or 1).
- **Output:** `ehh` list, decaying with distance.
- **Example code:**

```python
import numpy as np

from dnakit.variants import calculate_ehh_decay

rng = np.random.default_rng(8)
haplotypes = rng.integers(0, 2, size=(40, 12))
values = calculate_ehh_decay(haplotypes.tolist())["ehh"]
print([round(value, 6) for value in values[:4]])
```

- **Example output:**

```text
[0.454545, 0.242424, 0.121212, 0.030303]
```

## 34) Find Mendelian errors

- **Purpose:** Count genotypes in a pedigree that are inconsistent with Mendelian transmission.
- **Method:** Call scikit-allel's mendel_errors; the two parent columns are mother then father.
- **API:** `dnakit.variants.find_mendel_errors(parents: Sequence[Sequence[Sequence[int]]], progeny: Sequence[Sequence[Sequence[int]]])`
- **Input:** Parent genotypes (variants x 2 x 2) and progeny genotypes (variants x n x 2).
- **Output:** `errors` for each site and progeny, and `total_errors`.
- **Example code:**

```python
from dnakit.variants import find_mendel_errors

parents = [[[0, 0], [1, 1]], [[0, 1], [0, 1]], [[1, 1], [0, 0]]]
progeny = [[[0, 1]], [[1, 1]], [[0, 0]]]
result = find_mendel_errors(parents, progeny)
print(result["errors"], result["total_errors"])
```

- **Example output:**

```text
[[0], [0], [1]] 1
```

## 35) Phase a trio by transmission

- **Purpose:** Phase mother, father and child genotypes using Mendelian transmission.
- **Method:** Call scikit-allel's phase_by_transmission, resolving phase within windows of `window_size` sites where it can be determined.
- **API:** `dnakit.variants.phase_trio_by_transmission(genotypes: Sequence[Sequence[Sequence[int]]], window_size: int)`
- **Input:** Trio genotypes (variants x 3 x 2) in mother, father, child order; window size in sites.
- **Output:** `genotypes` (allele indices) and the `phased` flag matrix.
- **Example code:**

```python
from dnakit.variants import phase_trio_by_transmission

genotypes = [[[0, 1], [0, 0], [0, 1]], [[1, 1], [0, 1], [1, 1]], [[0, 0], [1, 0], [0, 0]], [[0, 1], [0, 1], [1, 1]]]
result = phase_trio_by_transmission(genotypes, 2)
print(result["genotypes"][0], result["phased"][0])
```

- **Example output:**

```text
[[0, 1], [0, 0], [1, 0]] [True, True, True]
```

## 36) Runs of homozygosity with an HMM

- **Purpose:** Call runs of homozygosity (ROH) from the genotypes of one diploid individual.
- **Method:** Call scikit-allel's roh_mhmm; a multinomial hidden Markov model separates ROH from non-ROH states using heterozygosity and transition probabilities.
- **API:** `dnakit.variants.call_roh_hmm(genotypes: Sequence[Sequence[int]], positions: Sequence[int], *, phet_roh: float = 0.001, transition: float = 1e-06, min_roh: int = 0, contig_size: int | None = None)`
- **Input:** One individual's genotypes (variants x 2), positions, heterozygosity and transition probabilities; optional `contig_size`.
- **Output:** `runs` (start, stop, length, is_marginal) and `froh`.
- **Example code:**

```python
import numpy as np

from dnakit.variants import call_roh_hmm

rng = np.random.default_rng(11)
calls = rng.choice(np.array([[0, 0], [0, 1]]), size=3000, p=[0.5, 0.5])
calls[1000:2000] = 0
positions = (np.arange(1, 3001) * 10).tolist()
result = call_roh_hmm(calls.tolist(), positions)
longest = max(result["runs"], key=lambda run: run["length"])
print(len(result["runs"]) > 0, longest["start"] <= 10100, longest["stop"] >= 19900)
```

- **Example output:**

```text
True True True
```

## 37) Regional nucleotide divergence dxy

- **Purpose:** Compute the average nucleotide divergence dxy between two populations over a region.
- **Method:** Sum the between-population mean differences over the sites in the region and divide by the number of bases (by default from the first to the last site).
- **API:** `dnakit.variants.calculate_sequence_divergence(positions: Sequence[int], first_counts: Sequence[Sequence[int]], second_counts: Sequence[Sequence[int]], *, start: int | None = None, stop: int | None = None)`
- **Input:** Positions and allele counts for two populations; optional region start and stop.
- **Output:** `dxy`, plus the region `start` and `stop`.
- **Example code:**

```python
from dnakit.variants import calculate_sequence_divergence

positions = [100, 200, 300, 400]
first = [[4, 0], [3, 1], [2, 2], [0, 4]]
second = [[0, 4], [1, 3], [2, 2], [4, 0]]
result = calculate_sequence_divergence(positions, first, second)
print(round(result["dxy"], 8), result["start"], result["stop"])
```

- **Example output:**

```text
0.01038206 100 400
```

## 38) Population branching statistic PBS

- **Purpose:** Compute PBS in sliding windows from three populations to locate lineage-specific differentiation.
- **Method:** Call scikit-allel's pbs, which combines three pairwise Fst values after a log transform; windows count sites.
- **API:** `dnakit.variants.calculate_population_branching_statistic(first_counts: Sequence[Sequence[int]], second_counts: Sequence[Sequence[int]], third_counts: Sequence[Sequence[int]], *, window_size: int, window_step: int | None = None, normed: bool = True)`
- **Input:** Allele counts for three populations, window size in sites, step and normalization flag.
- **Output:** `values` for each window.
- **Example code:**

```python
from dnakit.variants import calculate_population_branching_statistic

counts = [[3, 1], [2, 2], [4, 0], [1, 3], [0, 4], [2, 2], [3, 1], [1, 3], [2, 2], [4, 0]]
result = calculate_population_branching_statistic(counts, counts, counts, window_size=2)
print(result["values"])
```

- **Example output:**

```text
[0.0, 0.0, 0.0, 0.0, 0.0]
```

## 39) Normalize indel allele position

- **Purpose:** Shift an indel to the leftmost (or rightmost) equivalent position in a reference string.
- **Method:** Call hgvs's allele normalizer: trim shared suffixes and prefixes, then shift along the reference until no further shift is possible. Coordinates are 0-based and half-open.
- **API:** `dnakit.variants.normalize_indel_alleles(reference: str, start: int, stop: int, alleles: Sequence[str], *, left: bool = True, bound: int | None = None, ref_step: int = 1, shuffle: bool = True)`
- **Input:** Reference string, start and stop, two allele strings; optional direction `left`, bound `bound` and step `ref_step`.
- **Output:** Normalized `start`, `stop`, `alleles` and whether they changed (`changed`).
- **Example code:**

```python
from dnakit.variants import normalize_indel_alleles

left = normalize_indel_alleles("ACGTTTTA", 5, 6, ("T", ""))
right = normalize_indel_alleles("ACGTTTTA", 3, 4, ("T", ""), left=False, bound=8)
print(left["start"], left["stop"], right["start"], right["stop"])
```

- **Example output:**

```text
3 4 6 7
```

## 40) Trim the common prefix of allele strings

- **Purpose:** Remove the prefix shared by a list of allele strings and report how many bases were removed.
- **Method:** Call hgvs's prefix trimmer, which compares the lexicographic minimum and maximum of the strings.
- **API:** `dnakit.variants.trim_common_prefix_strings(strings: Sequence[str], *, min_length: int = 0)`
- **Input:** At least two allele strings; optional `min_length` to keep.
- **Output:** `trimmed_length` and the trimmed `strings`.
- **Example code:**

```python
from dnakit.variants import trim_common_prefix_strings

result = trim_common_prefix_strings(["ACGT", "ACTT"])
print(result["trimmed_length"], result["strings"])
```

- **Example output:**

```text
2 ['GT', 'TT']
```

## 41) Build a transcript CIGAR

- **Purpose:** Join per-exon CIGAR strings into one transcript CIGAR, writing introns as N.
- **Method:** Call hgvs's build_tx_cigar; on the minus strand each exon CIGAR is reversed, and intron length is the gap between neighbouring exons.
- **API:** `dnakit.variants.build_transcript_cigar(exons: Sequence[Mapping[str, Any]], strand: int)`
- **Input:** Exons in transcript order (cigar, alt_start_i, alt_end_i) and strand (1 or -1).
- **Output:** `cigar` string and `exon_count`.
- **Example code:**

```python
from dnakit.variants import build_transcript_cigar

exons = [
    {"cigar": "10=", "alt_start_i": 0, "alt_end_i": 10},
    {"cigar": "5=", "alt_start_i": 20, "alt_end_i": 25},
]
print(build_transcript_cigar(exons, 1)["cigar"])
```

- **Example output:**

```text
10=10N5=
```

## 42) Check accession fit for an HGVS type

- **Purpose:** Check whether an accession such as NM_ or NC_ fits an HGVS variant type (c, g, m, n, p or r).
- **Method:** Call hgvs's accession validator, which returns VALID, WARNING or ERROR from pattern matches.
- **API:** `dnakit.variants.validate_hgvs_type_accession(variant_type: str, accession: str)`
- **Input:** Variant type letter and accession string.
- **Output:** `level` (VALID, WARNING or ERROR) and `message`.
- **Example code:**

```python
from dnakit.variants import validate_hgvs_type_accession

print(validate_hgvs_type_accession("c", "NM_000088.3")["level"])
print(validate_hgvs_type_accession("g", "NM_000088.3")["level"])
```

- **Example output:**

```text
VALID
ERROR
```

## 43) Detect tandem repeat blocks

- **Purpose:** Find tandem repeat units in a sequence with their repeat counts.
- **Method:** Call hgvs's repeat block detector; scan left to right or right to left and report unit, count and block length.
- **API:** `dnakit.variants.detect_repeat_blocks(sequence: str, *, reverse: bool = False)`
- **Input:** Sequence (A/C/G/T/N only) and scan direction `reverse`.
- **Output:** `blocks`, each with repeat_unit, repeat_count, block_size and block.
- **Example code:**

```python
from dnakit.variants import detect_repeat_blocks

for block in detect_repeat_blocks("ACACACGTTTT")["blocks"]:
    print(block["repeat_unit"], block["repeat_count"], block["block_size"])
```

- **Example output:**

```text
AC 3 6
T 4 4
```

## 44) Assemble repeat blocks in HGVS notation

- **Purpose:** Format detected repeat blocks as HGVS-style repeat notation, such as AC[3]G[1]T[4].
- **Method:** Call hgvs's repeat string assembler; bases outside any repeat block are written as single-base counts.
- **API:** `dnakit.variants.assemble_repeat_blocks(sequence: str, blocks: Sequence[Mapping[str, Any]], *, reverse: bool = False)`
- **Input:** Original sequence, the block list from detect_repeat_blocks and direction `reverse`.
- **Output:** `formatted` string.
- **Example code:**

```python
from dnakit.variants import assemble_repeat_blocks, detect_repeat_blocks

sequence = "ACACACGTTTT"
blocks = detect_repeat_blocks(sequence)["blocks"]
print(assemble_repeat_blocks(sequence, blocks)["formatted"])
```

- **Example output:**

```text
AC[3]G[1]T[4]
```

## Limitations and notes

- varcode's `variant_overlaps_interval` ends a variant at `start + n_ref_bases`, one base past an inclusive end. The wrapper passes `n_ref_bases - 1` so the test follows the inclusive convention.
- scikit-allel's `heterozygosity_expected` uses ploidy as an exponent and gives the Hardy-Weinberg value only for diploids; the wrapper fixes ploidy at 2.
- Phase evidence for `enumerate_phase_hypotheses` and `partition_germline_by_phase` is supplied through `in_cis_evidence`; transcript-aware phasing needs pyensembl data and is not wrapped.
- The hgvs wrappers do not call seqrepo or UTA. Projection and alignment functions that need them (`full_house`, `variant_context_w_alignment`) are not wrapped.
- Missing genotype calls are −1 in input; missing statistics are None in output.
