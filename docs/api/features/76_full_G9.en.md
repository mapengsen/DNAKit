# Offline cutadapt and umi-tools sequence and UMI functions

These functions live in `dnakit.sequencing_qc` and wrap the Python interfaces of [cutadapt](https://github.com/marcelm/cutadapt) (5.x) and [umi-tools](https://github.com/CGATOxford/UMI-tools) (1.1.x). The cutadapt part computes neighbour sets, adapter k-mer search sets, poly-A and NextSeq trim positions, and adapter-specification and adapter-file parsing. The umi-tools part computes UMI distances and adjacency, barcode extraction and whitelists, CIGAR splice offsets, BAM tags and per-gene UMI counts. Every function works only on explicit sequences, files or count tables. None of them trims reads, runs the cutadapt or umi_tools command line, or deduplicates UMIs.

The dependencies are optional and are imported only when the matching function is called. If a package is missing, the function raises an error that includes the install command.

## Install optional dependencies

```bash
python -m pip install "dnakit[bio-cutadapt]"    # cutadapt and dnaio
python -m pip install "dnakit[bio-umi-tools]"   # umi-tools and its runtime imports
python -m pip install "dnakit[bio-reads]"       # pysam: BAM reading, and imported by umi-tools
```

Importing umi-tools also loads pysam, matplotlib, scipy and pybktree, and `detect_bam_tags()` also needs pandas. The `bio-umi-tools` extra currently lists only umi-tools and regex, so these packages must be installed separately (see Current scope).

Project links: [cutadapt](https://github.com/marcelm/cutadapt), [umi-tools](https://github.com/CGATOxford/UMI-tools), [pysam](https://github.com/pysam-developers/pysam).

Each example below shows output produced by actually running it. The examples import from the modules `dnakit.sequencing_qc` and `dnakit.sequencing_qc`; the package-level names in the API lines refer to the same functions.

## 1) `dnakit.sequencing_qc.hamming_neighbours()`

- **Purpose:** List every sequence whose Hamming distance from a DNA sequence is exactly (`exact=True`) or at most (default) `max_distance`. It is often used to build candidate sets for barcode correction.
- **Method:** With `exact=True` the function calls `cutadapt.align.hamming_sphere`. Otherwise it calls `cutadapt.align.hamming_environment` and returns distances 0 to `max_distance`. Each result gives its distance and number of matching bases. Results are sorted by distance and then by sequence. If the count exceeds `max_neighbours`, the function stops during enumeration. Only A, C, G and T are accepted.
- **API:** `dnakit.sequencing_qc.hamming_neighbours(sequence, max_distance, *, exact=False, max_neighbours=10000)`.
- **Input:** A DNA string; `max_distance` is an integer from 0 to 3 and not longer than the sequence; `exact` is a boolean; `max_neighbours` is an explicit limit.
- **Output:** A dictionary with `count` and `neighbours`. Each item in `neighbours` has `sequence`, `distance` and `matches`. `provenance` records the cutadapt version and the function used.
- **Example code:**

```python
from dnakit.sequencing_qc import hamming_neighbours

exact = hamming_neighbours("ACGT", 1, exact=True)
print(exact["count"], [row["sequence"] for row in exact["neighbours"]][:4])
within = hamming_neighbours("ACGT", 1)
print(within["count"], within["neighbours"][0])
```

- **Example output:**

```text
12 ['AAGT', 'ACAT', 'ACCT', 'ACGA']
13 {'sequence': 'ACGT', 'distance': 0, 'matches': 4}
```

## 2) `dnakit.sequencing_qc.edit_neighbours()`

- **Purpose:** List every sequence whose edit (Levenshtein) distance from a DNA sequence is at most `max_distance`. Insertions, deletions and substitutions all count.
- **Method:** Calls `cutadapt.align.edit_environment`, the Cython implementation. It gives each sequence its edit distance and the number of matches in the optimal alignment. A sequence keeps only its smallest distance. Results are sorted by distance and then by sequence.
- **API:** `dnakit.sequencing_qc.edit_neighbours(sequence, max_distance, *, max_neighbours=10000)`.
- **Input:** A DNA string; `max_distance` is an integer from 0 to 3; `max_neighbours` is an explicit limit.
- **Output:** `count`, `neighbours` (each item has `sequence`, `distance` and `matches`) and `provenance`.
- **Example code:**

```python
from dnakit.sequencing_qc import edit_neighbours

result = edit_neighbours("AC", 1)
print(result["count"], [row["sequence"] for row in result["neighbours"]][:5])
```

- **Example output:**

```text
19 ['AC', 'A', 'AA', 'AAC', 'ACA']
```

## 3) `dnakit.sequencing_qc.adapter_kmer_positions()`

- **Purpose:** Build cutadapt's k-mer search sets for one adapter. For each search offset it gives the k-mers that must occur in the read, which lets cutadapt locate the adapter quickly.
- **Method:** Calls `cutadapt.kmer_heuristic.create_positions_and_kmers`. `start` and `stop` are cutadapt offsets, where a negative value counts back from the read end. `stop` of `None` means the search extends to the read end. K-mers are sorted within each entry, and entries are sorted by `(start, stop)` so that results are deterministic.
- **API:** `dnakit.sequencing_qc.adapter_kmer_positions(adapter, *, min_overlap=3, max_error_rate=0.1, back_adapter=True, front_adapter=False, internal=True)`.
- **Input:** An adapter DNA string; `min_overlap` is a positive integer not longer than the adapter; `max_error_rate` is in [0, 1); the three booleans select the 3' adapter, the 5' adapter and internal matches.
- **Output:** `entries` (each with `start`, `stop`, `kmers`), `entry_count`, `kmer_count` and `provenance`.
- **Example code:**

```python
from dnakit.sequencing_qc import adapter_kmer_positions

result = adapter_kmer_positions("ACGTACGT", min_overlap=3)
for entry in result["entries"]:
    print(entry["start"], entry["stop"], entry["kmers"])
```

- **Example output:**

```text
-8 None ['ACGTA']
-4 None ['ACGT']
-3 None ['ACG']
0 None ['ACGTACGT']
```

## 4) `dnakit.sequencing_qc.adapter_kmer_table()`

- **Purpose:** Estimate how often the adapter's k-mer search sets are hit by random sequence. Use it to check whether an adapter's k-mers are too common.
- **Method:** Runs `create_positions_and_kmers` and then `cutadapt.kmer_heuristic.kmer_probability_analysis` with `default_length` as the read length. The upstream output is a tab-separated text table. This function parses each line into a record and extracts the overall hit probability. The raw text is kept in `table_text`.
- **API:** `dnakit.sequencing_qc.adapter_kmer_table(adapter, *, min_overlap=3, max_error_rate=0.1, back_adapter=True, front_adapter=False, internal=True, default_length=150)`.
- **Input:** The same inputs as `adapter_kmer_positions()`, plus `default_length` (an integer from 1 to 1000000).
- **Output:** `rows` (each with `kmer`, `start`, `stop`, `considered_sites`, `hit_chance_percent`), `profile_hit_chance_percent`, `table_text` and `provenance`.
- **Example code:**

```python
from dnakit.sequencing_qc import adapter_kmer_table

result = adapter_kmer_table("ACGTACGT")
for row in result["rows"]:
    print(row["kmer"], row["start"], row["considered_sites"], row["hit_chance_percent"])
print(result["profile_hit_chance_percent"])
```

- **Example output:**

```text
ACG -3 1 1.56
ACGT -4 1 0.39
ACGTA -8 4 0.39
ACGTACGT 0 143 0.22
2.54
```

## 5) `dnakit.sequencing_qc.sequence_kmer_chunks()`

- **Purpose:** Split a sequence into chunks of almost equal length and return the distinct chunk strings. cutadapt uses them as k-mer seeds.
- **Method:** Calls `cutadapt.kmer_heuristic.kmer_chunks`. The remainder of the sequence length divided by the chunk count decides how many chunks are one base longer. Results are sorted. IUPAC letters are accepted.
- **API:** `dnakit.sequencing_qc.sequence_kmer_chunks(sequence, chunks)`.
- **Input:** A sequence string; `chunks` is an integer from 1 to the sequence length.
- **Output:** `kmers` (a sorted list), `count`, `sequence`, `chunks` and `provenance`.
- **Example code:**

```python
from dnakit.sequencing_qc import sequence_kmer_chunks

result = sequence_kmer_chunks("ACGTACGTAC", 3)
print(result["kmers"])
```

- **Example output:**

```text
['ACG', 'ACGT', 'TAC']
```

## 6) `dnakit.sequencing_qc.poly_a_tail_start()`

- **Purpose:** Find the poly-A tail at the 3' end of a sequence. With `reverse_complement=True`, find the poly-T head at the 5' end instead.
- **Method:** Calls `cutadapt.qualtrim.poly_a_trim_index`. `index` is the start of the tail (or the end of the head), and `length` is the number of bases in the tail or head. cutadapt ignores tails shorter than 3 bases. With no tail, `index` equals the sequence length.
- **API:** `dnakit.sequencing_qc.poly_a_tail_start(sequence, *, reverse_complement=False)`.
- **Input:** A sequence string (IUPAC letters allowed) and the boolean `reverse_complement`.
- **Output:** `index`, `length`, `sequence_length`, `reverse_complement` and `provenance`.
- **Example code:**

```python
from dnakit.sequencing_qc import poly_a_tail_start

tail = poly_a_tail_start("ACGTAAAAAAAAA")
print(tail["index"], tail["length"])
head = poly_a_tail_start("TTTTTACGT", reverse_complement=True)
print(head["index"], head["length"])
```

- **Example output:**

```text
4 9
5 5
```

## 7) `dnakit.sequencing_qc.nextseq_trim_bounds()`

- **Purpose:** Compute the 3' trim position for NextSeq two-colour quality trimming. Runs of G at the 3' end, which carry no colour signal and often come from fragment ends, are removed.
- **Method:** Wraps the sequence and qualities in a `dnaio` record and calls `cutadapt.qualtrim.nextseq_trim_index`. The algorithm is the BWA running-sum rule used by `find_quality_trim_bounds()`, except that the quality of each G is taken as `cutoff - 1`. Only the 3' end is trimmed, so `index` is the retained length.
- **API:** `dnakit.sequencing_qc.nextseq_trim_bounds(sequence, qualities, *, cutoff=20, phred_offset=33)`.
- **Input:** A sequence; a Phred quality string of the same length; `cutoff` from 0 to 93; `phred_offset` of 33 or 64.
- **Output:** `index`, `retained_length`, `length`, `cutoff`, `phred_offset` and `provenance`.
- **Example code:**

```python
from dnakit.sequencing_qc import nextseq_trim_bounds

result = nextseq_trim_bounds("ACGTGGGGGG", "I" * 10, cutoff=20)
print(result["index"], result["retained_length"])
```

- **Example output:**

```text
4 4
```

## 8) `dnakit.sequencing_qc.expand_adapter_braces()`

- **Purpose:** Expand cutadapt's `x{n}` repeat notation, which makes adapters with repeated bases easier to write.
- **Method:** Calls `cutadapt.parser.expand_braces`, which replaces each `x{n}` with n copies of x. Before expansion the function estimates the length from the repeat counts and refuses inputs whose estimate exceeds `max_length`. Malformed notation, such as a missing number or unbalanced braces, raises `ConfigurationError`.
- **API:** `dnakit.sequencing_qc.expand_adapter_braces(sequence, *, max_length=1000000)`.
- **Input:** A string using `x{n}` notation; `max_length` is an explicit limit.
- **Output:** `expanded`, `input_length`, `expanded_length` and `provenance`.
- **Example code:**

```python
from dnakit.sequencing_qc import expand_adapter_braces

print(expand_adapter_braces("TGA{5}CT")["expanded"])
print(expand_adapter_braces("AC{3}G{2}")["expanded"])
```

- **Example output:**

```text
TGAAAAACT
ACCCGG
```

**Example file:** `adapters.fa`

```text
>ad1 description
ACGTACGT
>ad2
TTGGAACC
```

## 9) `dnakit.sequencing_qc.read_adapter_fasta()`

- **Purpose:** Read adapter names and sequences from a local FASTA file.
- **Method:** Calls `cutadapt.parser.read_adapters_fasta`. The name is the first word of the header, and an empty header gives `None`. Sequences are kept as written. Records are returned in file order. A file that is not FASTA, or that has more than `max_records` records, raises an error.
- **API:** `dnakit.sequencing_qc.read_adapter_fasta(source, *, max_records=10000)`.
- **Input:** A FASTA file path; `max_records` is an explicit limit.
- **Output:** `records` (each with `name` and `sequence`), `count`, `source` and `provenance`.
- **Example code:**

```python
from dnakit.sequencing_qc import read_adapter_fasta

result = read_adapter_fasta("adapters.fa")
for record in result["records"]:
    print(record["name"], record["sequence"])
```

- **Example output:**

```text
ad1 ACGTACGT
ad2 TTGGAACC
```

## 10) `dnakit.sequencing_qc.parse_adapter_spec()`

- **Purpose:** Parse one cutadapt adapter specification and return the adapter objects with their settings.
- **Method:** Calls `cutadapt.parser.make_adapters_from_one_specification`. A specification can be `ACGT`, `name=ACGT`, a linked pair `ACGT...TTTT` (not allowed with `anywhere`), or `file:<path>`, which reads every adapter in a local FASTA. After `file:`, parameters can be appended with semicolons, such as `;max_errors=0.2`. `adapter_type` is `front` (-g), `back` (-a) or `anywhere` (-b). A linked adapter reports its front and back sequences.
- **API:** `dnakit.sequencing_qc.parse_adapter_spec(spec, *, adapter_type="back", max_errors=0.1, min_overlap=3, indels=True, max_adapters=1000)`.
- **Input:** A specification string; `adapter_type`; `max_errors` in [0, 1); `min_overlap`; `indels` as a boolean; `max_adapters` as an explicit limit.
- **Output:** `adapters` (each with `name`, `adapter_class`, `sequence`, `max_error_rate`, `min_overlap`, `indels` and more; linked adapters have `front_sequence` and `back_sequence`), `count` and `provenance`.
- **Example code:**

```python
from dnakit.sequencing_qc import parse_adapter_spec

result = parse_adapter_spec("file:adapters.fa", adapter_type="back", max_errors=0.2)
for row in result["adapters"]:
    print(row["name"], row["adapter_class"], row["sequence"], row["max_error_rate"])
linked = parse_adapter_spec("ACGT...TTTT", adapter_type="back")
print(linked["adapters"][0]["front_sequence"], linked["adapters"][0]["back_sequence"])
```

- **Example output:**

```text
ad1 BackAdapter ACGTACGT 0.2
ad2 BackAdapter TTGGAACC 0.2
ACGT TTTT
```

## 11) `dnakit.sequencing_qc.average_umi_distance()`

- **Purpose:** Compute the average pairwise Hamming distance of a set of equal-length UMIs, to gauge UMI diversity.
- **Method:** Converts the UMIs to bytes and calls `umi_tools.umi_methods.get_average_umi_distance`, which needs bytes. Every pair is counted, so duplicate UMIs contribute distance 0. With a single UMI, upstream returns the sentinel -1 and this function reports `None`.
- **API:** `dnakit.sequencing_qc.average_umi_distance(umis, *, max_umis=2000)`.
- **Input:** A list of UMI strings (A, C, G, T, N) of equal length; `max_umis` is an explicit limit.
- **Output:** `umi_count`, `umi_length`, `pair_count`, `average_distance` and `provenance`.
- **Example code:**

```python
from dnakit.sequencing_qc import average_umi_distance

result = average_umi_distance(["AAAA", "AAAT", "CCCC"])
print(result["average_distance"], result["pair_count"])
```

- **Example output:**

```text
3.0 3
```

## 12) `dnakit.sequencing_qc.find_umi_neighbours()`

- **Purpose:** Find pairs of UMIs whose Hamming distance is at most `max_distance`, for screening candidates before UMI clustering.
- **Method:** Duplicates are removed and the UMIs are sorted first. `umi_tools.network.build_substr_idx` cuts each UMI into `max_distance + 1` pieces, so any two UMIs within the distance share at least one piece. `iter_nearest_neighbours` generates the candidate pairs, and the upstream `umi_tools.umi_methods.edit_distance` confirms each one. In each result `umi_a` sorts before `umi_b`, and the pairs are sorted.
- **API:** `dnakit.sequencing_qc.find_umi_neighbours(umis, *, max_distance=1, max_umis=2000)`.
- **Input:** A list of equal-length UMIs; `max_distance` is 0 to 3 and less than the UMI length.
- **Output:** `unique_umis`, `pairs` (each with `umi_a`, `umi_b`, `distance`), `pair_count` and `provenance`.
- **Example code:**

```python
from dnakit.sequencing_qc import find_umi_neighbours

result = find_umi_neighbours(["AAAA", "AAAT", "CCCC", "AAAG"], max_distance=1)
for pair in result["pairs"]:
    print(pair["umi_a"], pair["umi_b"], pair["distance"])
```

- **Example output:**

```text
AAAA AAAG 1
AAAA AAAT 1
AAAG AAAT 1
```

## 13) `dnakit.sequencing_qc.umi_connected_component()`

- **Purpose:** Return the connected component of an undirected UMI graph that contains a given node. It groups UMIs that are linked by mismatches.
- **Method:** Calls `umi_tools.network.breadth_first_search`. Neighbours that are missing as keys are added as nodes with no edges. The component is returned sorted.
- **API:** `dnakit.sequencing_qc.umi_connected_component(adjacency, node, *, max_nodes=100000)`.
- **Input:** `adjacency`, a mapping from each UMI to its neighbour list; `node`, the starting UMI; `max_nodes`, an explicit limit.
- **Output:** `component` (a sorted list), `size`, `node` and `provenance`.
- **Example code:**

```python
from dnakit.sequencing_qc import umi_connected_component

adjacency = {"A": ["B"], "B": ["A", "C"], "D": ["E"]}
result = umi_connected_component(adjacency, "A")
print(result["component"], result["size"])
```

- **Example output:**

```text
['A', 'B', 'C'] 3
```

## 14) `dnakit.sequencing_qc.mask_low_quality_umi()`

- **Purpose:** Replace UMI bases whose quality is below a threshold with N.
- **Method:** Calls `umi_tools.extract_methods.mask_umi`. The positions are found with `get_below_threshold`. `quality_encoding` is `phred33`, `phred64` or `solexa`. Only bases with a quality strictly below the threshold are masked.
- **API:** `dnakit.sequencing_qc.mask_low_quality_umi(umi, qualities, *, quality_encoding="phred33", threshold=10)`.
- **Input:** A UMI; a quality string of the same length; the encoding; a threshold from 0 to 93.
- **Output:** `masked_umi`, `masked_positions` (0-based), `masked_count` and `provenance`.
- **Example code:**

```python
from dnakit.sequencing_qc import mask_low_quality_umi

result = mask_low_quality_umi("ACGT", "II!I", threshold=10)
print(result["masked_umi"], result["masked_positions"])
```

- **Example output:**

```text
ACNT [2]
```

## 15) `dnakit.sequencing_qc.extract_barcode_positions()`

- **Purpose:** Remove UMI, cell-barcode and discarded bases from a read by 0-based position, and report the barcode sequences and qualities.
- **Method:** Calls `umi_tools.extract_methods.extractSeqAndQuals`. Upstream returns the retained qualities and the UMI and cell-barcode qualities, but not the barcode bases. This function reads the barcode sequences from the input at those positions, in read order. The three position sets must not overlap.
- **API:** `dnakit.sequencing_qc.extract_barcode_positions(sequence, qualities, *, umi_positions=(), cell_positions=(), discard_positions=(), retain_umi=False)`.
- **Input:** A read sequence; a Phred+33 quality string; three lists of 0-based positions; `retain_umi` as a boolean.
- **Output:** `read_sequence`, `read_qualities`, `umi_sequence`, `umi_qualities`, `cell_sequence`, `cell_qualities`, `positions` and `provenance`. Note: with `retain_umi=True` upstream copies the qualities of every retained base into `umi_qualities`; this function returns that value unchanged.
- **Example code:**

```python
from dnakit.sequencing_qc import extract_barcode_positions

result = extract_barcode_positions(
    "ACGTACGT", "IIIIHHHH",
    umi_positions=[0, 1], cell_positions=[2], discard_positions=[7],
)
print(result["read_sequence"], result["read_qualities"])
print(result["umi_sequence"], result["umi_qualities"], result["cell_sequence"])
```

- **Example output:**

```text
TACG IHHH
AC II G
```

## 16) `dnakit.sequencing_qc.tag_read_name_with_umi()`

- **Purpose:** Append a UMI, and optionally a cell barcode, after the first word of a FASTQ header.
- **Method:** Calls `umi_tools.extract_methods.addBarcodesToIdentifier`. Upstream splits the header on spaces and appends `separator + [cell barcode + separator] + UMI` to the first word. The rest of the header is unchanged.
- **API:** `dnakit.sequencing_qc.tag_read_name_with_umi(read_name, umi, *, cell_barcode=None, separator="_")`.
- **Input:** A header string; a UMI; an optional cell barcode; a separator of 1 to 4 non-space characters.
- **Output:** `tagged_name`, `umi`, `cell_barcode`, `separator` and `provenance`.
- **Example code:**

```python
from dnakit.sequencing_qc import tag_read_name_with_umi

result = tag_read_name_with_umi("read1 extra words", "ACGT", cell_barcode="TTGG")
print(result["tagged_name"])
```

- **Example output:**

```text
read1_TTGG_ACGT extra words
```

## 17) `dnakit.sequencing_qc.cigar_splice_offset()`

- **Purpose:** Return the offset, from the read start, of the first splice junction (N) or trailing soft clip (S) in a CIGAR string. It supports splice-junction checks in RNA-seq.
- **Method:** Parses the CIGAR into `(operation code, length)` pairs and calls `umi_tools.sam_methods.find_splice`. The offset counts only operations that consume the reference (M, D, =, X), and skips insertions, hard clips and padding. A leading soft clip is added to the offset but is not a splice. `reverse=True` scans from the 3' end, as umi-tools does for reverse-strand reads. `splice_offset` is `None` when there is no splice.
- **API:** `dnakit.sequencing_qc.cigar_splice_offset(cigar, *, reverse=False)`.
- **Input:** A CIGAR string such as `10M100N20M`; the boolean `reverse`.
- **Output:** `splice_offset`, `has_splice`, `reverse`, `cigar` and `provenance`.
- **Example code:**

```python
from dnakit.sequencing_qc import cigar_splice_offset

print(cigar_splice_offset("10M100N20M")["splice_offset"])
print(cigar_splice_offset("30M")["splice_offset"])
print(cigar_splice_offset("10M100N20M", reverse=True)["splice_offset"])
```

- **Example output:**

```text
10
None
20
```

## 18) `dnakit.sequencing_qc.parse_umi_read_name()`

- **Purpose:** Read the UMI, and optionally the cell barcode, from the end of a read name.
- **Method:** Calls `umi_tools.sam_methods.get_umi_read_string` (UMI only) or `get_cell_umi_read_string` (cell barcode, then UMI). Both split the name on the separator and use the last field, and for the cell barcode the field before it. Too few fields raise an error.
- **API:** `dnakit.sequencing_qc.parse_umi_read_name(read_name, *, separator="_", cell_barcode=False)`.
- **Input:** A read name; the separator; whether to read the cell barcode too.
- **Output:** `umi`, `cell_barcode` (`None` when not requested), `mode` and `provenance`.
- **Example code:**

```python
from dnakit.sequencing_qc import parse_umi_read_name

result = parse_umi_read_name("read1_CELL1_UMI1", cell_barcode=True)
print(result["cell_barcode"], result["umi"])
print(parse_umi_read_name("read1_UMI1")["umi"])
```

- **Example output:**

```text
CELL1 UMI1
UMI1
```

## 19) `dnakit.sequencing_qc.barcode_indel_matches()`

- **Purpose:** Find whitelist barcodes that differ from a barcode by a few substitutions, insertions or deletions. It is used for barcode correction.
- **Method:** Calls `umi_tools.whitelist_methods.checkError`, which uses fuzzy regex matching. Upstream skips a whitelist entry identical to the barcode and stops after two matches. So `matches` holds at most two barcodes, and `ambiguous` is true when at least two were found. `exact_in_whitelist` is true when the barcode itself is in the whitelist; that entry is not reported as a match.
- **API:** `dnakit.sequencing_qc.barcode_indel_matches(barcode, whitelist, *, max_errors=1)`.
- **Input:** A barcode (A, C, G, T); a whitelist list; `max_errors` is an integer from 0 to 3.
- **Output:** `matches`, `match_count`, `ambiguous`, `exact_in_whitelist`, `whitelist_size` and `provenance`.
- **Example code:**

```python
from dnakit.sequencing_qc import barcode_indel_matches

result = barcode_indel_matches("ACGT", ["ACGA", "ACGTT", "TTTT"])
print(result["matches"], result["ambiguous"])
```

- **Example output:**

```text
['ACGA', 'ACGTT'] True
```

**Example file:** `whitelist.tsv`

```text
# barcode<TAB>corrections
AAAA	CCCC,CCCG
GGGG	TTTT
```

## 20) `dnakit.sequencing_qc.read_barcode_whitelist()`

- **Purpose:** Read a user-supplied cell-barcode whitelist from a TSV file, with an optional error-correction map.
- **Method:** Calls `umi_tools.whitelist_methods.getUserDefinedBarcodes`. The first column of each line is a whitelist barcode. Lines starting with `#` are comments, and blank lines raise an error. `error_correction` is `none` (the whitelist only), `file` (the second column lists comma-separated barcodes that correct to the first column), or `derived` (each base of a whitelist barcode is replaced by another character from A, C, G, T and N, and `threshold`, 1 or 2, is the maximum number of replaced bases). With `threshold=2`, a barcode within range of two whitelist barcodes is left out of the map.
- **API:** `dnakit.sequencing_qc.read_barcode_whitelist(source, *, error_correction="none", threshold=1, max_barcodes=100000)`.
- **Input:** A TSV path; `error_correction`; `threshold` (1 or 2); `max_barcodes`, an explicit limit.
- **Output:** `barcodes`, `barcode_count`, `error_map` (`None` without correction), `error_correction`, `threshold` and `provenance`. `error_map` maps each corrected barcode to its whitelist barcode.
- **Example code:**

```python
from dnakit.sequencing_qc import read_barcode_whitelist

result = read_barcode_whitelist("whitelist.tsv", error_correction="file")
print(result["barcodes"])
print(result["error_map"])
```

- **Example output:**

```text
['AAAA', 'GGGG']
{'CCCC': 'AAAA', 'CCCG': 'AAAA', 'TTTT': 'GGGG'}
```

## 21) `dnakit.sequencing_qc.call_cell_whitelist()`

- **Purpose:** Decide which barcodes are real cell barcodes from their read counts, using a knee-point method.
- **Method:** Calls `umi_tools.whitelist_methods.getCellWhitelist`. `knee_method="distance"` takes the point of the cumulative curve farthest from the line joining its first and last points (`getKneeEstimateDistance`). `"density"` uses a Gaussian density (`getKneeEstimateDensity`). `expect_cells` applies only to density, and `cell_number` fixes how many barcodes are accepted. With `error_correct_threshold` above 0, the result includes, for each accepted barcode, the non-whitelist barcodes within that Hamming distance (`error_map`). Upstream needs a `collections.Counter`, which this function builds.
- **API:** `dnakit.sequencing_qc.call_cell_whitelist(cell_barcode_counts, *, knee_method="distance", expect_cells=None, cell_number=None, error_correct_threshold=0, max_barcodes=100000)`.
- **Input:** A mapping from barcode to count (positive integers); `knee_method`; optional `expect_cells` and `cell_number`; `error_correct_threshold` (0 to 2).
- **Output:** `whitelist`, `cell_count`, `input_barcodes`, `knee_method`, `error_map` (or `None`) and `provenance`.
- **Example code:**

```python
from dnakit.sequencing_qc import call_cell_whitelist

counts = {"AAAA": 100, "AAAT": 90, "CCCC": 10, "GGGG": 9, "TTTT": 1, "ACGT": 50}
result = call_cell_whitelist(counts, knee_method="distance")
print(result["whitelist"], result["cell_count"])
```

- **Example output:**

```text
['AAAA', 'AAAT', 'ACGT'] 3
```

**Example file:** `reads.sam`

```text
@HD	VN:1.6	SO:unsorted
@SQ	SN:chr1	LN:100
r1	0	chr1	1	60	4M	*	0	0	ACGT	IIII	NH:i:1
r2	0	chr1	9	60	4M	*	0	0	ACGT	IIII	NH:i:1	XT:A:U
```

## 22) `dnakit.sequencing_qc.detect_bam_tags()`

- **Purpose:** Check whether the first mapped reads of a SAM or BAM file carry the NH, X0 and XT tags.
- **Method:** Calls `umi_tools.dedup.detect_bam_features`, which reads the file through pysam. Only mapped reads are examined. A tag is reported `True` only when every examined read carries it. Upstream examines `n_entries + 1` records.
- **API:** `dnakit.sequencing_qc.detect_bam_tags(source, *, n_entries=1000)`.
- **Input:** A `.bam` or `.sam` path; `n_entries` is an integer from 1 to 1000000.
- **Output:** `tags` (boolean values for NH, X0 and XT), `n_entries`, `source` and `provenance`.
- **Example code:**

```python
from dnakit.sequencing_qc import detect_bam_tags

result = detect_bam_tags("reads.sam")
print(result["tags"])
```

- **Example output:**

```text
{'NH': True, 'X0': False, 'XT': False}
```

**Example file:** `assignments.tsv`

```text
r1_CELL1_UMI1	GENE_A
r2_CELL1_UMI1	GENE_A
r3_CELL1_UMI2	GENE_A
r4_CELL2_UMI1	GENE_B
r5_CELL1_UMI1	GENE_B
```

## 23) `dnakit.sequencing_qc.count_umis_per_gene()`

- **Purpose:** Count reads per gene, and the number of distinct UMIs (or distinct cell-barcode and UMI pairs) per gene.
- **Method:** The input is a headerless two-column TSV of read name and assigned gene. Calls `umi_tools.sam_methods.get_gene_count_tab`, which takes the UMI from the read name (and the cell barcode too when `cell_barcode=True`). Upstream expects rows for one gene to be contiguous. This function merges the groups itself, so the input order does not matter. Upstream yields a `None` gene for empty input, and this function ignores it.
- **API:** `dnakit.sequencing_qc.count_umis_per_gene(source, *, cell_barcode=False, separator="_", max_rows=1000000)`.
- **Input:** A TSV file with two non-empty columns per line; the separator in read names; `max_rows`, an explicit limit.
- **Output:** `genes` (sorted by gene; each with `gene`, `reads`, `unique_umis`), `gene_count`, `row_count`, `cell_barcode` and `provenance`.
- **Example code:**

```python
from dnakit.sequencing_qc import count_umis_per_gene

result = count_umis_per_gene("assignments.tsv", cell_barcode=True)
for gene in result["genes"]:
    print(gene["gene"], gene["reads"], gene["unique_umis"])
```

- **Example output:**

```text
GENE_A 3 2
GENE_B 2 2
```

## How these functions are verified

The test file is `others/tests/unit/test_full_G9.py`, with 30 tests. They check the following:

- Every wrapped function has at least one test. Where possible the test compares with a direct call to upstream, for example `cutadapt.align.hamming_sphere`, `cutadapt.qualtrim.poly_a_trim_index`, `umi_tools.network.breadth_first_search` and `umi_tools.whitelist_methods.getCellWhitelist`. The others are checked against hand calculations. For example, 12 sequences are exactly one Hamming step from `ACGT`, and 19 sequences are within edit distance 1 of `AC`. The 19 is also confirmed by exhaustive search with upstream `edit_distance`.
- For UMI pair search, random UMIs are compared with an exhaustive search that uses upstream Hamming distance, which confirms that every pair within the distance is found.
- File-based functions (adapter FASTA, whitelist, BAM, gene count table) use small files in `tmp_path` and do not depend on repository data.
- Invalid input (non-DNA letters, out-of-range or overlapping positions, values over a limit, malformed adapter braces) raises `ConfigurationError`.
- With `cutadapt` or `umi_tools` hidden through `sys.modules`, the tests check that `BackendUnavailableError` names `bio-cutadapt` or `bio-umi-tools`.
- Grouped results (k-mer search sets, UMI graphs, whitelist maps, gene counts) are identical when `PYTHONHASHSEED` is 1 and when it is 2.
- No function accesses the network, so there are no network functions to mock.

## Current scope

- This batch wraps 23 functions: 10 from cutadapt and 13 from umi-tools. The other entries of the two upstream packages are classified in the table below (the counts come from this batch's ledger records). Every function accepts the keyword-only argument `verbose=True` for progress output, which is omitted from the API lines.

| Status | Meaning | cutadapt | umi-tools |
| --- | --- | ---: | ---: |
| `wrapped` (Wrapped) | Wrapped in this batch | 10 | 13 |
| `already-wrapped` (Already wrapped) | Already wrapped in `dnakit.sequencing_qc` | 4 | 2 |
| `duplicate-of` (Duplicate) | Another dnakit function already gives the same computation | 3 | 5 |
| `family-variant` (Family variant) | Covered by a parameter or alias of a wrapped function | 13 | 7 |
| `cli` (Command-line entry) | Command-line or script entry point; not re-wrapped | 22 | 13 |
| `internal` (Internal helper) | Logging, temporary files, resource limits or internal steps | 27 | 16 |
| `other` (Other) | Needs cutadapt report objects, pysam records or pandas data | 11 | 11 |

- Relation to existing wrappers: `find_quality_trim_bounds()` and `calculate_expected_errors()` (cutadapt `qualtrim`), and `extract_read_barcodes()` and `map_barcodes_to_whitelist()` (umi-tools), already exist in `dnakit.sequencing_qc`, so this batch does not wrap them again. `find_adapter_match()` and `trim_fastq_adapters()` wrap cutadapt classes such as `BackAdapter` and `AdapterCutter`, which have no function entries in the inventory.
- pyslow5: the inventory has no function entries. pyslow5 1.5.0 is a single extension module with no submodules. Its public interface is the class `pyslow5.Open`, with methods such as `seq_reads`, `get_read_list` and `write_record`. The existing `read_slow5_signal()`, `summarize_slow5_signals()` and `write_slow5_signals()` live in `dnakit.sequencing_qc`. **Caution: pyslow5 1.5.0 can crash the interpreter on corrupt SLOW5/BLOW5 files.** The existing wrappers check the file's magic bytes before opening it.
- Dependencies: umi-tools modules load pysam, matplotlib, scipy and pybktree when they are imported, and `umi_tools.dedup` also loads pandas, but the `bio-umi-tools` extra lists only umi-tools and regex. The functions in this batch check these packages before importing the upstream modules, and a missing one reports the `bio-umi-tools` hint, so that extra needs to be completed.
- Upstream problems (recorded, not fixed):
  1. `umi_tools.umi_methods.edit_distance` raises `AttributeError` for inputs of unequal length, because upstream uses `np.Inf`, which NumPy 2 removed. This batch calls it only on UMIs of one length.
  2. With `retain_umi=True`, `extractSeqAndQuals` copies the qualities of every retained base into `umi_qualities`, not only those of UMI bases. This function returns that value unchanged.
  3. The whitelist and knee functions require `collections.Counter`; a plain dictionary raises `AttributeError`. This function builds the Counter itself.
  4. `detect_bam_features` opens its pysam file without closing it, so the handle stays open until garbage collection.
  5. `get_gene_count_tab` expects rows for each gene to be contiguous, and empty input yields a `None` gene. This function merges the groups and ignores `None`.
  6. `checkError` stops after two matches, so `matches` is not a complete list. This function reports `ambiguous` in its result.
