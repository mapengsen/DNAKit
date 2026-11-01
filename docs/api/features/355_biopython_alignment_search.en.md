# Biopython alignment formats, MAF index, HMMER and Infernal reports and alignment analysis

This group of functions in `dnakit.alignment.biopython_p18_alignment_functions` wraps [Biopython](https://biopython.org/) (`Bio.Align`, `Bio.AlignIO.MafIO`, `Bio.SearchIO`, `Bio.Nexus`) for the alignment operations that the earlier Biopython wrappers did not reach. They read and write the alignment formats A2M, bigBed, bigMaf, bigPsl, chain, Clustal, EMBOSS, Exonerate, MAF, Mauve, MSF, NEXUS, PHYLIP, PSL, Stockholm and tabular output; query a MAF file by reference region; read HMMER and Infernal reports and index, filter and sort search results; run the FOGSAA and codon-aware pairwise aligners; and analyse alignments (row sorting, substitution counts, column frequencies, identity counts, coordinate maps, merging, mapping) and DNA matrices (constant and gap-only columns, site states, subsetting, gap insertion).

Alignments are passed as a mapping of sequence name to gapped sequence (`-` for a gap) and returned the same way, so every result is JSON-serialisable. Writers never overwrite: the output path must be a new file. Biopython is an optional dependency (`pip install "dnakit[bio-formats]"`) imported only when a function is called. The examples import from the module, `dnakit.alignment.biopython_p18_alignment_functions`.

Not included on purpose: `QueryResult.hit_map`, `hsp_map` and `Hit.map` apply an arbitrary Python function to every hit or HSP. A function object cannot be passed through the JSON-typed interface, so the operation is skipped; `filter_search_results()` and `sort_search_results()` cover the declarative cases.

## Example data

Three declared inputs are used: a MAF file with two alignment blocks on the reference `hg.chr1`, a reduced HMMER `--domtblout` table with two queries, and an Infernal `cmscan --tblout` table with three queries (the HMMER tabular formats cannot be indexed by Biopython, so the indexing example uses the Infernal file).

**Example file:** `demo.maf`

```text
##maf version=1 scoring=none
a score=100
s hg.chr1 100 10 + 1000 ACGTACGTAC
s mm.chr2 200  9 + 2000 ACGTA-GTAC

a score=80
s hg.chr1 120  8 + 1000 GGTTAACC
s mm.chr2 400  8 + 2000 GGTTAACC
```

**Example file:** `hits.domtab`

```text
#                                                                            --- full sequence --- -------------- this domain -------------   hmm coord   ali coord   env coord
# target name        accession   tlen query name           accession   qlen   E-value  score  bias   #  of  c-Evalue  i-Evalue  score  bias  from    to  from    to  from    to  acc description of target
#------------------- ---------- ----- -------------------- ---------- ----- --------- ------ ----- --- --- --------- --------- ------ ----- ----- ----- ----- ----- ----- ----- ---- ---------------------
Globin               PF00042.17   110 gi|4885477|ref|NP_005359.1| -            154     1e-22   80.5   0.3   1   1   1.1e-26   1.6e-22   79.8   0.3     1   109     7   112     7   113 0.97 Globin
Xpo1                 PF08389.7    148 gi|22748937|ref|NP_065801.1| -           1204   8.5e-34  116.6   7.8   1   2   2.3e-37   1.2e-33  116.1   4.9     2   148   110   271   109   271 0.98 Exportin 1-like protein
Xpo1                 PF08389.7    148 gi|22748937|ref|NP_065801.1| -           1204   8.5e-34  116.6   7.8   2   2      0.52   2.6e+03   -1.8   0.0   112   139   499   525   496   529 0.86 Exportin 1-like protein
IBN_N                PF03810.14    77 gi|22748937|ref|NP_065801.1| -           1204    0.0044   16.9   0.0   1   2   7.3e-06     0.036   13.9   0.0     4    75    36    98    33   100 0.87 Importin-beta N-terminal domain
IBN_N                PF03810.14    77 gi|22748937|ref|NP_065801.1| -           1204    0.0044   16.9   0.0   2   2         2     1e+04   -3.5   0.0    57    74   168   185   165   186 0.84 Importin-beta N-terminal domain
Rac1                 PF09632.5    297 gi|22748937|ref|NP_065801.1| -           1204     0.095   11.7   0.3   1   2   0.00026       1.3    8.0   0.3    94   171    89   169    82   254 0.75 Rac1-binding domain
Rac1                 PF09632.5    297 gi|22748937|ref|NP_065801.1| -           1204     0.095   11.7   0.3   2   2     0.041     2e+02    0.8   0.0   150   181   310   341   296   354 0.82 Rac1-binding domain
#
# [ok]
```

**Example file:** `cm.tbl`

```text
#target name         accession query name              accession mdl mdl from   mdl to seq from   seq to strand trunc pass   gc  bias  score   E-value inc description of target
#------------------- --------- ----------------------- --------- --- -------- -------- -------- -------- ------ ----- ---- ---- ----- ------ --------- --- ---------------------
U2                   RF00004   ENA|BK006935|BK006935.2 -          cm        1      193    52929    53083      +    no    1 0.44   0.0   13.5      0.91 ?   U2 spliceosomal RNA
U2                   RF00004   ENA|BK006935|BK006935.2 -          cm        1      193   196571   196389      -    no    1 0.33   5.3   12.8       1.3 ?   U2 spliceosomal RNA
U2                   RF00004   ENA|BK006936|BK006936.2 -          cm        1      193   681858   681747      -    no    1 0.33   0.1   98.7   1.2e-20 !   U2 spliceosomal RNA
5S_rRNA              RF00001   ENA|BK006937|BK006937.2 -          cm        1      119      761      644      -    no    1 0.41   0.3   14.1       2.4 ?   5S ribosomal RNA
U2                   RF00004   ENA|BK006937|BK006937.2 -          cm        1      193   229986   229885      -    no    1 0.32   0.1   11.1       4.7 ?   U2 spliceosomal RNA
#
# Program:         cmscan
# Version:         1.1.5 (Sep 2023)
# Pipeline mode:   SCAN
# Query file:      BK006935_6_7.fasta
# Target file:     IRES_5S_U2.cm
# Option settings: cmscan --tblout IRES_5S_U2_Yeast-cmscan.tbl IRES_5S_U2.cm BK006935_6_7.fasta 
# Current dir:     /analysis/BioPython/Testing
# Date:            Thu Sep 12 05:06:09 2024
# [ok]
```

## 1) `dnakit.alignment.biopython_p18_alignment_functions.read_alignment_records()`

- **Purpose:** Read alignments of a format that the other alignment readers do not cover.
- **Method:** Call `Bio.Align.parse()` for `a2m`, `bigmaf` (bigMaf, a binary bigBed file), `bigpsl` (bigPsl), `chain` (UCSC chain), `emboss` (EMBOSS needle/water/matcher output), `exonerate` (text, cigar and vulgar output are detected), `mauve` (XMFA), `msf`, `nexus`, `stockholm` or `tabular` (BLAST or FASTA tabular output written with `-outfmt 7`). Block coordinates are 0-based and half-open, one row per sequence.
- **API:** `dnakit.alignment.biopython_p18_alignment_functions.read_alignment_records(path, *, file_format='nexus', max_alignments=10000)`.
- **Input:** A local file path, `file_format` and `max_alignments` (default 10000).
- **Output:** `alignment_count` and `alignments`; each has `sequence_ids`, `sequence_lengths`, `aligned_columns`, `score`, `coordinates`, the gapped `aligned_rows` when the file stores sequences, `annotations` and `column_annotations`.
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.alignment.biopython_p18_alignment_functions import (
    read_alignment_records,
    write_aligned_sequences,
)

path = Path(tempfile.mkdtemp()) / "demo.nex"
write_aligned_sequences({"s1": "ACGT-ACGTA", "s2": "AC-TTACGAA"}, path, file_format="nexus")
result = read_alignment_records(path, file_format="nexus")
alignment = result["alignments"][0]
print(result["alignment_count"], alignment["sequence_ids"], alignment["aligned_columns"])
print(alignment["coordinates"])
print(alignment["aligned_rows"])
```

- **Example output:**

```text
1 ['s1', 's2'] 10
[[0, 2, 3, 4, 4, 9], [0, 2, 2, 3, 4, 9]]
{'s1': 'ACGT-ACGTA', 's2': 'AC-TTACGAA'}
```

## 2) `dnakit.alignment.biopython_p18_alignment_functions.write_aligned_sequences()`

- **Purpose:** Write one multiple alignment (sequence name to gapped sequence) to a new file.
- **Method:** Build a `Bio.Align.Alignment` from the gapped rows (`Alignment.parse_printed_alignment()`) and write it with `Bio.Align.write()` in `a2m`, `clustal`, `fasta`, `maf`, `mauve` (XMFA), `nexus`, `phylip` or `stockholm` format. `molecule_type` (`DNA` or `RNA`) is used by the NEXUS writer. PHYLIP names are padded or limited by Biopython's writer.
- **API:** `dnakit.alignment.biopython_p18_alignment_functions.write_aligned_sequences(aligned_sequences, output_path, *, file_format='clustal', molecule_type='DNA')`.
- **Input:** A mapping of name to gapped sequence (equal lengths, names without spaces), a new output path, `file_format` and `molecule_type`.
- **Output:** `output_path`, `file_format`, `sequence_count`, `aligned_columns` and `alignment_count`.
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.alignment.biopython_p18_alignment_functions import write_aligned_sequences

out = Path(tempfile.mkdtemp()) / "demo.phy"
rows = {"s1": "ACGT-ACGTA", "s2": "AC-TTACGAA", "s3": "ACGTTACG-A"}
result = write_aligned_sequences(rows, out, file_format="phylip")
print(result["sequence_count"], result["aligned_columns"], result["alignment_count"])
print(out.read_text(encoding="utf-8"))
```

- **Example output:**

```text
3 10 1
3 10
s1        ACGT-ACGTA
s2        AC-TTACGAA
s3        ACGTTACG-A
```

## 3) `dnakit.alignment.biopython_p18_alignment_functions.convert_alignment_file()`

- **Purpose:** Convert an alignment file with `Bio.Align` to A2M, bigBed, bigMaf, bigPsl, chain, Clustal, Exonerate, MAF, Mauve, PHYLIP, PSL or Stockholm.
- **Method:** Parse the source with `Bio.Align.parse()` (sources: a2m, bigbed, bigmaf, bigpsl, chain, clustal, emboss, exonerate, fasta, maf, mauve, msf, nexus, phylip, psl, stockholm, tabular) and write it with `Bio.Align.write()`. `exonerate_variant` is `vulgar` or `cigar` for the Exonerate target. The text targets that print sequences (a2m, clustal, maf, mauve, phylip, stockholm) need a source that stores sequences; the coordinate targets (bigBed, bigMaf, bigPsl, chain, Exonerate, PSL) need a pairwise source that Biopython can write in that format. The attributes that Biopython's writers read from alignments of another format (target sizes, annotations, column states, Mauve identifiers) are supplied, and alignments are sorted by target name and start for the bigBed family. If Biopython refuses a combination, the function stops with its message and removes the partial file.
- **API:** `dnakit.alignment.biopython_p18_alignment_functions.convert_alignment_file(source, output_path, *, source_format='psl', target_format='chain', exonerate_variant='vulgar', max_alignments=10000)`.
- **Input:** A source path, a new output path, `source_format`, `target_format`, `exonerate_variant` and `max_alignments`.
- **Output:** `alignment_count`, `bytes_written`, `source_format` and `target_format`.
- **Example code:**

```python
import tempfile
from pathlib import Path

from Bio import Align
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from dnakit.alignment.biopython_p18_alignment_functions import convert_alignment_file

aligner = Align.PairwiseAligner(
    mode="global", match_score=2, mismatch_score=-1, open_gap_score=-2, extend_gap_score=-1
)
target = SeqRecord(Seq("ACGTACGTTGCA"), id="chr1")
query = SeqRecord(Seq("ACGTCGTTTGCA"), id="read1")
work = Path(tempfile.mkdtemp())
Align.write(aligner.align(target, query)[0], work / "demo.psl", "psl")
result = convert_alignment_file(work / "demo.psl", work / "demo.chain", source_format="psl", target_format="chain")
print(result["alignment_count"], result["target_format"])
print((work / "demo.chain").read_text(encoding="utf-8"))
```

- **Example output:**

```text
1 chain
chain 0 chr1 12 + 0 12 read1 12 + 0 12
4	1	0
2	0	1
5
```

## 4) `dnakit.alignment.biopython_p18_alignment_functions.search_maf_index()`

- **Purpose:** Find the alignment blocks of a MAF file that overlap regions of the reference sequence.
- **Method:** Build a `Bio.AlignIO.MafIO.MafIndex` over `target_seqname` (for example `mm9.chr10`) and call `search()`. The index lives in the SQLite file `index_path` and is reused when it exists; without `index_path` a temporary index is used and removed. `starts` are 0-based and `ends` half-open, as in the UCSC convention.
- **API:** `dnakit.alignment.biopython_p18_alignment_functions.search_maf_index(maf_path, target_seqname, starts, ends, *, index_path=None, max_blocks=1000)`.
- **Input:** A MAF path, `target_seqname`, `starts`, `ends`, optional `index_path` and `max_blocks`.
- **Output:** `indexed_blocks`, `block_count` and `blocks`; each block lists its `rows` (`id`, gapped `sequence`, MAF `annotations` such as start, size, strand and srcSize).
- **Example code:**

```python
from dnakit.alignment.biopython_p18_alignment_functions import search_maf_index

result = search_maf_index("demo.maf", "hg.chr1", [102], [108])
print(result["indexed_blocks"], result["block_count"])
for row in result["blocks"][0]["rows"]:
    print(row["id"], row["sequence"], row["annotations"]["start"])
```

- **Example output:**

```text
2 1
hg.chr1 ACGTACGTAC 100
mm.chr2 ACGTA-GTAC 200
```

## 5) `dnakit.alignment.biopython_p18_alignment_functions.splice_maf_index()`

- **Purpose:** Extract the in-silico spliced alignment of several reference segments, for example exons.
- **Method:** Call `MafIndex.get_spliced()` on the same index as `search_maf_index()`: the segments `starts`/`ends` (0-based, half-open) on `target_seqname` are cut out of the MAF alignment and joined into one multiple alignment with a row per species. `strand` is 1 or -1; -1 reverse-complements the result.
- **API:** `dnakit.alignment.biopython_p18_alignment_functions.splice_maf_index(maf_path, target_seqname, starts, ends, *, strand=1, index_path=None)`.
- **Input:** A MAF path, `target_seqname`, `starts`, `ends`, `strand` and optional `index_path`.
- **Output:** `alignment_length`, `segment_count`, `strand` and `rows` (sequence name to spliced gapped sequence).
- **Example code:**

```python
from dnakit.alignment.biopython_p18_alignment_functions import splice_maf_index

forward = splice_maf_index("demo.maf", "hg.chr1", [100, 120], [104, 124])
print(forward["alignment_length"], forward["rows"])
reverse = splice_maf_index("demo.maf", "hg.chr1", [100, 120], [104, 124], strand=-1)
print(reverse["rows"])
```

- **Example output:**

```text
8 {'hg.chr1': 'ACGTGGTT', 'mm.chr2': 'ACGTGGTT'}
{'hg.chr1': 'AACCACGT', 'mm.chr2': 'AACCACGT'}
```

## 6) `dnakit.alignment.biopython_p18_alignment_functions.read_search_results()`

- **Purpose:** Read a search-program report (HMMER, Infernal and other formats) with `Bio.SearchIO`.
- **Method:** Call `Bio.SearchIO.parse()` for `hmmer2-text`, `hmmer3-text`, `hmmer3-tab` (`--tblout`), `hmmscan3-domtab`, `hmmsearch3-domtab`, `phmmer3-domtab` (`--domtblout`), `infernal-tab`, `infernal-text`, `blast-tab`, `blast-xml`, `blat-psl`, `exonerate-text`, `exonerate-cigar`, `exonerate-vulgar` or `fasta-m10`. `comments=True` reads BLAST tabular files written with `-outfmt 7`. HSP ranges are 0-based and half-open, as Bio.SearchIO stores them.
- **API:** `dnakit.alignment.biopython_p18_alignment_functions.read_search_results(path, *, file_format='hmmer3-text', comments=False, max_queries=10000)`.
- **Input:** A local report path, `file_format`, `comments` and `max_queries`.
- **Output:** `query_count`, `hit_count` and `queries`; each query has its program, target, version, `hits` (id, description, E-value and bit score where the format gives them) and their `hsps` (ranges, scores, strands and the format's own attributes).
- **Example code:**

```python
from dnakit.alignment.biopython_p18_alignment_functions import read_search_results

result = read_search_results("hits.domtab", file_format="hmmscan3-domtab")
print(result["query_count"], result["hit_count"])
for query in result["queries"]:
    for hit in query["hits"]:
        first = hit["hsps"][0]
        print(hit["id"], hit["hsp_count"], first["evalue"], first["query_range"])
```

- **Example output:**

```text
2 4
Globin 1 1.6e-22 [6, 112]
Xpo1 2 1.2e-33 [109, 271]
IBN_N 2 0.036 [35, 98]
Rac1 2 1.3 [88, 169]
```

## 7) `dnakit.alignment.biopython_p18_alignment_functions.index_search_results()`

- **Purpose:** Index a search-output file for random access by query id and return chosen queries.
- **Method:** Without `index_path` an in-memory offset index is built (`Bio.SearchIO.index()`); with it the index is stored in a SQLite file (`Bio.SearchIO.index_db()`), created when missing and reused otherwise. The queries named in `query_ids` are read through the index. Biopython 1.88 cannot index the HMMER tabular formats (`hmmer3-tab` and the three `domtab` formats); they are refused with a clear message.
- **API:** `dnakit.alignment.biopython_p18_alignment_functions.index_search_results(path, *, file_format='hmmer3-text', index_path=None, query_ids=None, comments=False)`.
- **Input:** A report path, `file_format`, optional `index_path`, `query_ids` and `comments`.
- **Output:** `query_count`, the sorted `query_ids`, `index_path` and the requested `queries` in full.
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.alignment.biopython_p18_alignment_functions import index_search_results

database = Path(tempfile.mkdtemp()) / "hits.sqlite"
result = index_search_results(
    "cm.tbl", file_format="infernal-tab", index_path=database,
    query_ids=["ENA|BK006936|BK006936.2"],
)
print(result["query_count"], result["query_ids"])
print(result["queries"][0]["id"], result["queries"][0]["hit_count"], database.exists())
```

- **Example output:**

```text
3 ['ENA|BK006935|BK006935.2', 'ENA|BK006936|BK006936.2', 'ENA|BK006937|BK006937.2']
ENA|BK006936|BK006936.2 1 True
```

## 8) `dnakit.alignment.biopython_p18_alignment_functions.filter_search_results()`

- **Purpose:** Filter the hits and HSPs of a search report.
- **Method:** HSPs are kept by `QueryResult.hsp_filter()` when their E-value is at most `max_evalue`, their bit score at least `min_bitscore` and their alignment span at least `min_alignment_span`; hits that lose all HSPs disappear. Then `QueryResult.hit_filter()` keeps hits with at least `min_hsps_per_hit` HSPs and, if `hit_ids` is given, whose id is in that list. Unset criteria are ignored; at least one must be given.
- **API:** `dnakit.alignment.biopython_p18_alignment_functions.filter_search_results(path, *, file_format='hmmer3-text', max_evalue=None, min_bitscore=None, min_alignment_span=None, min_hsps_per_hit=None, hit_ids=None, comments=False)`.
- **Input:** A report path, `file_format` and the criteria above.
- **Output:** `hits_before`, `hits_after` and the remaining `queries`.
- **Example code:**

```python
from dnakit.alignment.biopython_p18_alignment_functions import filter_search_results

strict = filter_search_results("hits.domtab", file_format="hmmscan3-domtab", max_evalue=1e-5)
print(strict["hits_before"], strict["hits_after"])
print([hit["id"] for query in strict["queries"] for hit in query["hits"]])
multi = filter_search_results("hits.domtab", file_format="hmmscan3-domtab", min_hsps_per_hit=2)
print(multi["hits_after"])
```

- **Example output:**

```text
4 2
['Globin', 'Xpo1']
3
```

## 9) `dnakit.alignment.biopython_p18_alignment_functions.sort_search_results()`

- **Purpose:** Sort the hits of every query, and the HSPs of every hit, of a search report.
- **Method:** `QueryResult.sort()` and `Hit.sort()` with a named key: `evalue` (ascending; hits without an E-value last), `bitscore`, `id`, `length` (hit sequence length where the format stores it) or `hsp_count`. `reverse` flips the order. HSPs are sorted by the same key when `sort_hsps` is true and the key applies to HSPs (evalue, bitscore).
- **API:** `dnakit.alignment.biopython_p18_alignment_functions.sort_search_results(path, *, file_format='hmmer3-text', sort_by='evalue', reverse=False, sort_hsps=True, comments=False)`.
- **Input:** A report path, `file_format`, `sort_by`, `reverse` and `sort_hsps`.
- **Output:** `sort_by`, `reverse` and the sorted `queries`.
- **Example code:**

```python
from dnakit.alignment.biopython_p18_alignment_functions import sort_search_results

result = sort_search_results("hits.domtab", file_format="hmmscan3-domtab", sort_by="bitscore", reverse=True)
query = result["queries"][1]
print([(hit["id"], hit["bitscore"]) for hit in query["hits"]])
print([hsp["bitscore"] for hsp in query["hits"][0]["hsps"]])
```

- **Example output:**

```text
[('Xpo1', 116.6), ('IBN_N', 16.9), ('Rac1', 11.7)]
[116.1, -1.8]
```

## 10) `dnakit.alignment.biopython_p18_alignment_functions.fogsaa_global_alignment()`

- **Purpose:** Align two sequences globally with the FOGSAA algorithm.
- **Method:** Use `Bio.Align.PairwiseAligner(mode='fogsaa')`, a branch-and-bound search for the optimal global alignment. The scores default to Biopython's (match 1, mismatch 0, no gap penalty); pass explicit scores for a useful DNA alignment.
- **API:** `dnakit.alignment.biopython_p18_alignment_functions.fogsaa_global_alignment(sequence_a, sequence_b, *, match_score=1.0, mismatch_score=0.0, open_gap_score=0.0, extend_gap_score=0.0, max_alignments=1)`.
- **Input:** `sequence_a`, `sequence_b`, `match_score`, `mismatch_score`, `open_gap_score`, `extend_gap_score` and `max_alignments`.
- **Output:** `score`, `alignment_count`, `parameters` and `alignments` (gapped rows named `sequence_a` and `sequence_b`, block coordinates).
- **Example code:**

```python
from dnakit.alignment.biopython_p18_alignment_functions import fogsaa_global_alignment

result = fogsaa_global_alignment(
    "ACGTACGTTGCA", "ACGTCGTTTGCA",
    match_score=2, mismatch_score=-1, open_gap_score=-2, extend_gap_score=-1,
)
print(result["score"], result["alignment_count"])
rows = result["alignments"][0]["aligned_rows"]
print(rows["sequence_a"])
print(rows["sequence_b"])
```

- **Example output:**

```text
18.0 1
ACGTACGTT-GCA
ACGT-CGTTTGCA
```

## 11) `dnakit.alignment.biopython_p18_alignment_functions.align_nucleotide_to_protein()`

- **Purpose:** Align a coding nucleotide sequence to its protein with the codon-aware aligner.
- **Method:** `Bio.Align.CodonAligner` pairs each amino acid with its codon and tolerates frame shifts and introns of the nucleotide sequence; `anchor_length` is the seed length used to find the next codon after a gap. DNA or RNA is accepted (RNA is read as DNA). `codon_table` is an NCBI genetic-code id. The score is the number of aligned codons.
- **API:** `dnakit.alignment.biopython_p18_alignment_functions.align_nucleotide_to_protein(protein, nucleotides, *, codon_table=1, anchor_length=10, max_alignments=10)`.
- **Input:** `protein`, `nucleotides`, `codon_table`, `anchor_length` and `max_alignments`.
- **Output:** `score`, `alignment_count` and `alignments` with block coordinates (protein in residues, nucleotide in bases).
- **Example code:**

```python
from dnakit.alignment.biopython_p18_alignment_functions import align_nucleotide_to_protein

plain = align_nucleotide_to_protein("MAIVMG", "ATGGCCATTGTAATGGGC")
print(plain["score"], plain["alignments"][0]["coordinates"])
shifted = align_nucleotide_to_protein(
    "MLGSNESRVCPKLSPS", "ATGCTGGGCTCGAACGAGTCCGTGTATGCCCTAAGCTGAGCCCGTCG"
)
print(shifted["score"], shifted["alignments"][0]["coordinates"])
```

- **Example output:**

```text
6.0 [[0, 6], [0, 18]]
13.0 [[0, 7, 7, 16], [0, 21, 20, 47]]
```

## 12) `dnakit.alignment.biopython_p18_alignment_functions.sort_alignment_rows()`

- **Purpose:** Sort the rows of an alignment.
- **Method:** Build an `Alignment` and call `Alignment.sort()` with a named key: `id` (the Biopython default), `sequence` (the ungapped letters), `gc_fraction` (of the ungapped sequence), `ungapped_length` or `gap_count`. Ties keep the input order.
- **API:** `dnakit.alignment.biopython_p18_alignment_functions.sort_alignment_rows(aligned_sequences, *, sort_by='id', reverse=False)`.
- **Input:** A mapping of name to gapped sequence, `sort_by` and `reverse`.
- **Output:** `order` (names in the new order) and `aligned_sequences` in that order.
- **Example code:**

```python
from dnakit.alignment.biopython_p18_alignment_functions import sort_alignment_rows

rows = {"r1": "ACGT-ACGTA", "r2": "AC-TTACGAA", "r3": "ACGTTACG-A"}
print(sort_alignment_rows(rows, sort_by="id", reverse=True)["order"])
print(sort_alignment_rows(rows, sort_by="gc_fraction", reverse=True)["order"])
print(sort_alignment_rows(rows, sort_by="gap_count", reverse=True)["order"])
```

- **Example output:**

```text
['r3', 'r2', 'r1']
['r1', 'r3', 'r2']
['r1', 'r2', 'r3']
```

## 13) `dnakit.alignment.biopython_p18_alignment_functions.alignment_substitution_counts()`

- **Purpose:** Count letter substitutions between the rows of an alignment.
- **Method:** Read `Alignment.substitutions`: entry (row letter, column letter) is the number of columns where the first letter is aligned to the second one (symmetrised for alignments of several rows, as Biopython does). Gap columns are not counted.
- **API:** `dnakit.alignment.biopython_p18_alignment_functions.alignment_substitution_counts(aligned_sequences)`.
- **Input:** A mapping of name to gapped sequence (at least two rows).
- **Output:** `alphabet` and the square `matrix`.
- **Example code:**

```python
from dnakit.alignment.biopython_p18_alignment_functions import alignment_substitution_counts

result = alignment_substitution_counts({"a": "ACGTTACGA", "b": "ACGATACGA", "c": "ACCTTACGA"})
print(result["alphabet"])
for letter, row in zip(result["alphabet"], result["matrix"]):
    print(letter, row)
```

- **Example output:**

```text
ACGT
A [9.0, 0.0, 0.0, 1.0]
C [0.0, 6.0, 0.0, 0.0]
G [0.0, 2.0, 4.0, 0.0]
T [1.0, 0.0, 0.0, 4.0]
```

## 14) `dnakit.alignment.biopython_p18_alignment_functions.alignment_column_frequencies()`

- **Purpose:** Count each letter per alignment column.
- **Method:** Read `Alignment.frequencies`: for each letter (and `-` for a gap) the number of rows that carry it in each column.
- **API:** `dnakit.alignment.biopython_p18_alignment_functions.alignment_column_frequencies(aligned_sequences)`.
- **Input:** A mapping of name to gapped sequence (at least two rows).
- **Output:** `column_count`, `row_count` and `frequencies` (letter to a list of counts, one per column).
- **Example code:**

```python
from dnakit.alignment.biopython_p18_alignment_functions import alignment_column_frequencies

rows = {"r1": "ACGT-ACGTA", "r2": "AC-TTACGAA", "r3": "ACGTTACG-A"}
result = alignment_column_frequencies(rows)
print(result["row_count"], result["column_count"])
for letter, counts in result["frequencies"].items():
    print(letter, counts)
```

- **Example output:**

```text
3 10
- [0.0, 0.0, 1.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0]
A [3.0, 0.0, 0.0, 0.0, 0.0, 3.0, 0.0, 0.0, 1.0, 3.0]
C [0.0, 3.0, 0.0, 0.0, 0.0, 0.0, 3.0, 0.0, 0.0, 0.0]
G [0.0, 0.0, 2.0, 0.0, 0.0, 0.0, 0.0, 3.0, 0.0, 0.0]
T [0.0, 0.0, 0.0, 3.0, 2.0, 0.0, 0.0, 0.0, 1.0, 0.0]
```

## 15) `dnakit.alignment.biopython_p18_alignment_functions.alignment_identity_counts()`

- **Purpose:** Count identities, mismatches and gaps of an alignment.
- **Method:** Call `Alignment.counts()`. Returned are the aligned columns, identities, mismatches, gaps and their split into left, internal and right gaps and into insertions and deletions. `wildcard` (a single character such as `N`) excludes that letter from the counts; `scoring_matrix` (a bundled substitution matrix name such as `NUC.4.4`) adds the `substitution_score` and the number of `positives`.
- **API:** `dnakit.alignment.biopython_p18_alignment_functions.alignment_identity_counts(aligned_sequences, *, wildcard=None, scoring_matrix=None)`.
- **Input:** A mapping of name to gapped sequence, `wildcard` or `scoring_matrix`.
- **Output:** `row_count`, `column_count` and `counts`.
- **Example code:**

```python
from dnakit.alignment.biopython_p18_alignment_functions import alignment_identity_counts

pair = {"a": "ACGT-ACGTA", "b": "AC-TTACGAA"}
counts = alignment_identity_counts(pair)["counts"]
print(counts["aligned"], counts["identities"], counts["mismatches"], counts["gaps"], counts["internal_gaps"])
scored = alignment_identity_counts(pair, scoring_matrix="NUC.4.4")["counts"]
print(scored["positives"], scored["substitution_score"])
```

- **Example output:**

```text
8 7 1 2 2
7 31.0
```

## 16) `dnakit.alignment.biopython_p18_alignment_functions.alignment_coordinate_maps()`

- **Purpose:** Map between alignment columns and sequence positions.
- **Method:** This covers `Alignment.aligned`, `Alignment.indices` and `Alignment.inverse_indices`. `indices` has, for each row and column, the 0-based position in the ungapped sequence or -1 at a gap; `inverse_indices` has, for each ungapped position, its column; `aligned_blocks` (two rows only) lists the aligned block intervals per sequence and is null for more rows.
- **API:** `dnakit.alignment.biopython_p18_alignment_functions.alignment_coordinate_maps(aligned_sequences)`.
- **Input:** A mapping of name to gapped sequence (at least two rows).
- **Output:** `indices`, `inverse_indices`, `aligned_blocks` and `column_count`.
- **Example code:**

```python
from dnakit.alignment.biopython_p18_alignment_functions import alignment_coordinate_maps

result = alignment_coordinate_maps({"a": "AC-GT", "b": "ACAG-"})
print(result["indices"])
print(result["inverse_indices"])
print(result["aligned_blocks"])
```

- **Example output:**

```text
{'a': [0, 1, -1, 2, 3], 'b': [0, 1, 2, 3, -1]}
{'a': [0, 1, 3, 4], 'b': [0, 1, 2, 3]}
[[[0, 2], [2, 3]], [[0, 2], [3, 4]]]
```

## 17) `dnakit.alignment.biopython_p18_alignment_functions.parse_printed_alignment()`

- **Purpose:** Infer the ungapped sequences and block coordinates from printed alignment lines.
- **Method:** Call `Alignment.parse_printed_alignment()`: one gapped line per sequence, all with the same number of columns and `-` for a gap. Returned are the ungapped sequences and the coordinate matrix (rows are sequences, columns are block boundaries, 0-based), which together define an `Alignment`.
- **API:** `dnakit.alignment.biopython_p18_alignment_functions.parse_printed_alignment(aligned_lines)`.
- **Input:** A list of gapped lines.
- **Output:** `sequences` and `coordinates`.
- **Example code:**

```python
from dnakit.alignment.biopython_p18_alignment_functions import parse_printed_alignment

result = parse_printed_alignment(["TAGGCATACGTG", "AACG--TACGT-", "-ACGCATACTTG"])
print(result["sequences"])
for row in result["coordinates"]:
    print(row)
```

- **Example output:**

```text
['TAGGCATACGTG', 'AACGTACGT', 'ACGCATACTTG']
[0, 1, 4, 6, 11, 12]
[0, 1, 4, 4, 9, 9]
[0, 0, 3, 5, 10, 11]
```

## 18) `dnakit.alignment.biopython_p18_alignment_functions.merge_alignments_by_reference()`

- **Purpose:** Combine pairwise alignments that share one reference sequence into a multiple alignment.
- **Method:** Call `Alignment.from_alignments_with_same_reference()`. Each query maps to `{"target": gapped reference, "query": gapped query}`; after removing gaps every reference must be the same sequence. The reference comes first in the result and columns are added where a query inserts letters.
- **API:** `dnakit.alignment.biopython_p18_alignment_functions.merge_alignments_by_reference(reference_alignments)`.
- **Input:** A mapping of query name to a `target`/`query` pair of gapped rows.
- **Output:** `row_count`, `column_count` and `aligned_sequences` (the reference is named `reference`).
- **Example code:**

```python
from dnakit.alignment.biopython_p18_alignment_functions import merge_alignments_by_reference

result = merge_alignments_by_reference(
    {
        "q1": {"target": "ACGT", "query": "AC-T"},
        "q2": {"target": "ACGT-", "query": "A-GTT"},
    }
)
print(result["row_count"], result["column_count"])
for name, row in result["aligned_sequences"].items():
    print(name, row)
```

- **Example output:**

```text
3 5
reference ACGT-
q1 AC-T-
q2 A-GTT
```

## 19) `dnakit.alignment.biopython_p18_alignment_functions.map_alignment_through()`

- **Purpose:** Map (lift over) an alignment through one or more alignments of its query sequence.
- **Method:** `alignment` is `{"target": gapped row, "query": gapped row}`, for example a chromosome against a transcript. `through` lists alignments of the same form whose target is the first alignment's query (ungapped), for example reads aligned to the transcript. One alignment uses `Alignment.map()`; several use `Alignment.mapall()`. The result aligns the mapped sequences to the first target.
- **API:** `dnakit.alignment.biopython_p18_alignment_functions.map_alignment_through(alignment, through)`.
- **Input:** The first alignment and a list of alignments to map through.
- **Output:** `method` (`map` or `mapall`), `coordinates` and the gapped `aligned_rows`.
- **Example code:**

```python
from dnakit.alignment.biopython_p18_alignment_functions import map_alignment_through

chromosome = "AAAAAAAACCCCCCCAAAAAAAAAAAGGGGGGAAAAAAAA"
transcript_alignment = {"target": chromosome, "query": "--------CCCCCCC-----------GGGGGG--------"}
read_alignment = {"target": "CCCCCCCGGGGGG", "query": "---CCCCGGGG--"}
result = map_alignment_through(transcript_alignment, [read_alignment])
print(result["method"], result["coordinates"])
print(result["aligned_rows"]["target"])
print(result["aligned_rows"]["query"])
```

- **Example output:**

```text
map [[11, 15, 26, 30], [0, 4, 4, 8]]
CCCCAAAAAAAAAAAGGGG
CCCC-----------GGGG
```

## 20) `dnakit.alignment.biopython_p18_alignment_functions.find_alignment_columns()`

- **Purpose:** List the constant or the gap-only columns of a DNA alignment.
- **Method:** `kind="constant"` calls `Nexus.constant()`: the columns where all (undeleted) sequences agree, treating `?` as compatible and ambiguity codes as sets of bases. `kind="gap_only"` calls `Nexus.gaponly()`: the columns that hold only gaps, counting `?` as a gap when `include_missing` is true. `delete_taxa` leaves sequences out and `exclude_sites` leaves columns out of the constant search. Columns are 0-based.
- **API:** `dnakit.alignment.biopython_p18_alignment_functions.find_alignment_columns(aligned_sequences, *, kind='constant', delete_taxa=None, exclude_sites=None, include_missing=False)`.
- **Input:** A mapping of name to gapped DNA, `kind`, `delete_taxa`, `exclude_sites` and `include_missing`.
- **Output:** `kind`, `columns`, `found` and `column_count`.
- **Example code:**

```python
from dnakit.alignment.biopython_p18_alignment_functions import find_alignment_columns

dna = {"t1": "ACGT-ACG", "t2": "ACGT-ATG", "t3": "ACGA-ACG"}
print(find_alignment_columns(dna)["columns"])
print(find_alignment_columns(dna, kind="gap_only")["columns"])
print(find_alignment_columns(dna, delete_taxa=["t3"])["columns"])
```

- **Example output:**

```text
[0, 1, 2, 4, 5, 7]
[4]
[0, 1, 2, 3, 4, 5, 7]
```

## 21) `dnakit.alignment.biopython_p18_alignment_functions.site_character_states()`

- **Purpose:** Summarise the character states of alignment columns.
- **Method:** Call `Nexus.cstatus()` for each column in `sites` (all columns when null; 0-based): the sorted list of bases seen in the undeleted sequences, with ambiguity codes expanded. `narrow=True` follows PAUP (`A C ?` gives `A C`; all `?` gives `?`); `narrow=False` expands `?` and gaps to all states.
- **API:** `dnakit.alignment.biopython_p18_alignment_functions.site_character_states(aligned_sequences, *, sites=None, delete_taxa=None, narrow=True)`.
- **Input:** A mapping of name to gapped DNA, `sites`, `delete_taxa` and `narrow`.
- **Output:** `sites`, a list of `{site, states}`, and `column_count`.
- **Example code:**

```python
from dnakit.alignment.biopython_p18_alignment_functions import site_character_states

dna = {"t1": "ACGT-ACG", "t2": "ACGT-ATG", "t3": "ACGA-ACG"}
for item in site_character_states(dna, sites=[3, 4, 6])["sites"]:
    print(item["site"], item["states"])
print(site_character_states(dna, sites=[4], narrow=False)["sites"][0]["states"])
```

- **Example output:**

```text
3 ['A', 'T']
4 ['?']
6 ['C', 'T']
['-', 'A', 'C', 'G', 'T']
```

## 22) `dnakit.alignment.biopython_p18_alignment_functions.subset_alignment()`

- **Purpose:** Delete sequences and exclude columns from a DNA alignment.
- **Method:** Call `Nexus.crop_matrix()`. `delete_taxa` names the sequences to drop and `exclude_sites` the 0-based columns to remove. At least one of them must be given and at least one sequence must remain.
- **API:** `dnakit.alignment.biopython_p18_alignment_functions.subset_alignment(aligned_sequences, *, delete_taxa=None, exclude_sites=None)`.
- **Input:** A mapping of name to gapped DNA, `delete_taxa` and `exclude_sites`.
- **Output:** `row_count`, `column_count` and `aligned_sequences`.
- **Example code:**

```python
from dnakit.alignment.biopython_p18_alignment_functions import subset_alignment

dna = {"t1": "ACGT-ACG", "t2": "ACGT-ATG", "t3": "ACGA-ACG"}
result = subset_alignment(dna, delete_taxa=["t2"], exclude_sites=[0, 1])
print(result["row_count"], result["column_count"])
print(result["aligned_sequences"])
```

- **Example output:**

```text
2 6
{'t1': 'GT-ACG', 't3': 'GA-ACG'}
```

## 23) `dnakit.alignment.biopython_p18_alignment_functions.insert_alignment_gaps()`

- **Purpose:** Insert gap columns into a DNA alignment.
- **Method:** Call `Nexus.insert_gap()`: `count` columns of `-` are added to every sequence before the 0-based column `position`. When the insertion point lies inside a run of gaps, the new gaps are placed at the right end of the run, or at the left end with `left_greedy=True`.
- **API:** `dnakit.alignment.biopython_p18_alignment_functions.insert_alignment_gaps(aligned_sequences, position, *, count=1, left_greedy=False)`.
- **Input:** A mapping of name to gapped DNA, `position`, `count` and `left_greedy`.
- **Output:** `column_count` and `aligned_sequences`.
- **Example code:**

```python
from dnakit.alignment.biopython_p18_alignment_functions import insert_alignment_gaps

dna = {"t1": "ACGT-ACG", "t2": "ACGT-ATG", "t3": "ACGA-ACG"}
result = insert_alignment_gaps(dna, 2, count=2)
print(result["column_count"])
for name, row in result["aligned_sequences"].items():
    print(name, row)
```

- **Example output:**

```text
10
t1 AC--GT-ACG
t2 AC--GT-ATG
t3 AC--GA-ACG
```
