# scikit-bio file formats, sequence queries and alignment helpers

These functions in `dnakit.sequence.skbio_p20_io_functions` are adapters over [scikit-bio](https://scikit.bio/): readers and writers for QSEQ, Clustal, PHYLIP, Newick, distance-matrix, ordination, QIIME 2 metadata, BIOM and NCBI taxdump files, sequence queries (masking, positions, motifs, reverse complements) and alignment helpers (consensus, conservation, gap frequencies, CIGAR paths). scikit-bio is an optional dependency (`dnakit[bio-skibio]`; BIOM needs `biom-format`) imported only when a function runs. Sequences must be upper-case IUPAC text as scikit-bio requires; positions are 0-based and intervals half-open. Tested with scikit-bio 0.7.4.

## 1) `dnakit.sequence.skbio_p20_io_functions.read_skbio_file()`

- **Purpose:** Read a bioinformatics file with scikit-bio and return its content as plain lists and mappings.
- **Method:** `file_format` selects the scikit-bio reader: `qseq` (Illumina QSEQ reads), `clustal` and `phylip` (alignments), `newick` (tree), `lsmat`, `phylip_dm` and `binary_dm` (distance matrices), `ordination`, `sample_metadata` (QIIME 2 TSV), `biom` and `taxdump` (NCBI `*.dmp`, with `taxdump_scheme` naming the columns). One function covers these eleven reader rows of the work list. Lists are cut to `max_items`.
- **API:** `dnakit.sequence.skbio_p20_io_functions.read_skbio_file(path, file_format, *, sequence_type='DNA', qseq_variant='illumina1.8', taxdump_scheme='nodes', max_items=1000)`.
- **Input:** A local file path and the format name; `sequence_type` (DNA or RNA) types alignment residues, `qseq_variant` is the quality encoding.
- **Output:** `format`, plus format-specific keys such as `alignment`, `newick`, `ids` and `matrix`, `records`, `ordination`, and provenance.
- **Example code:**

```python
from pathlib import Path
Path("qseq_reads.txt").write_text("M1\t1\t1\t1\t1\t1\t1\t1\tACG.\tIIII\t1\n")
Path("aln.clustal").write_text("CLUSTAL W (1.83)\n\ns1 ACGT-A\ns2 ACGTTA\n\n")
Path("aln.phy").write_text(" 2 6\ns1        ACGT-A\ns2        ACGTTA\n")
Path("tree.nwk").write_text("((a:1,b:2)90:1,c:3);")
Path("dm.phylip").write_text("2\na 0 1\nb 1 0\n")
from dnakit.sequence.skbio_p20_io_functions import read_skbio_file
reads = read_skbio_file("qseq_reads.txt", "qseq", verbose=False)
print(reads["records"][0]["sequence"], reads["records"][0]["quality"])
for name, fmt in (("aln.clustal", "clustal"), ("aln.phy", "phylip")):
    print(fmt, read_skbio_file(name, fmt, verbose=False)["alignment"])
tree = read_skbio_file("tree.nwk", "newick", verbose=False)
print(tree["newick"], tree["tip_count"])
dm = read_skbio_file("dm.phylip", "phylip_dm", verbose=False)
print(dm["ids"], dm["matrix"])
```

- **Example output:**

```text
ACG. [40, 40, 40, 40]
clustal {'s1': 'ACGT-A', 's2': 'ACGTTA'}
phylip {'s1': 'ACGT-A', 's2': 'ACGTTA'}
((a:1.0,b:2.0)90:1.0,c:3.0); 3
['a', 'b'] [[0.0, 1.0], [1.0, 0.0]]
```

## 2) `dnakit.sequence.skbio_p20_io_functions.write_skbio_distance_matrix()`

- **Purpose:** Write a symmetric distance matrix as lsmat text or binary (HDF5).
- **Method:** Builds a scikit-bio `DistanceMatrix` and writes it with the `lsmat` or `binary_dm` writer; `file_format` picks the format.
- **API:** `dnakit.sequence.skbio_p20_io_functions.write_skbio_distance_matrix(distance_matrix, output_path, *, ids=None, file_format='lsmat')`.
- **Input:** A square list matrix, optional `ids`, an output path.
- **Output:** The written path, ids and size.
- **Example code:**

```python
from dnakit.sequence.skbio_p20_io_functions import read_skbio_file, write_skbio_distance_matrix
matrix = [[0, 1, 2], [1, 0, 1.5], [2, 1.5, 0]]
for fmt, name in (("lsmat", "d.lsmat"), ("binary_dm", "d.h5")):
    out = write_skbio_distance_matrix(matrix, name, ids=["a", "b", "c"], file_format=fmt, verbose=False)
    back = read_skbio_file(name, fmt, verbose=False)
    print(out["format"], out["size"], back["ids"], back["matrix"][2])
```

- **Example output:**

```text
lsmat 3 ['a', 'b', 'c'] [2.0, 1.5, 0.0]
binary_dm 3 ['a', 'b', 'c'] [2.0, 1.5, 0.0]
```

## 3) `dnakit.sequence.skbio_p20_io_functions.write_skbio_ordination()`

- **Purpose:** Write ordination results (sample coordinates and eigenvalues) in scikit-bio's ordination text format.
- **Method:** Builds an `OrdinationResults` object from the coordinates and eigenvalues and writes it with the `ordination` writer.
- **API:** `dnakit.sequence.skbio_p20_io_functions.write_skbio_ordination(sample_coordinates, eigenvalues, output_path, *, sample_ids=None, proportion_explained=None, method='PCoA', long_method_name='Principal Coordinate Analysis')`.
- **Input:** Sample-by-axis coordinates, one eigenvalue per axis, optional sample ids and proportion explained, an output path.
- **Output:** The written path, sample count and axis count.
- **Example code:**

```python
from dnakit.sequence.skbio_p20_io_functions import read_skbio_file, write_skbio_ordination
out = write_skbio_ordination([[0.1, 0.2], [0.3, -0.2], [-0.4, 0.0]], [2.0, 1.0], "pcoa.ord",
                             sample_ids=["a", "b", "c"], proportion_explained=[0.6, 0.4], verbose=False)
print(out["sample_count"], out["axes"])
back = read_skbio_file("pcoa.ord", "ordination", verbose=False)["ordination"]
print(back["samples"]["index"], back["proportion_explained"])
```

- **Example output:**

```text
3 2
['a', 'b', 'c'] {'0': 0.6, '1': 0.4}
```

## 4) `dnakit.sequence.skbio_p20_io_functions.write_skbio_sample_metadata()`

- **Purpose:** Write a QIIME 2 sample metadata TSV with a type row.
- **Method:** Builds a `SampleMetadata` object from the columns (numbers become numeric, strings categorical) and calls `SampleMetadata.save`.
- **API:** `dnakit.sequence.skbio_p20_io_functions.write_skbio_sample_metadata(columns, sample_ids, output_path)`.
- **Input:** A mapping of column name to one value per sample, the sample ids and an output path.
- **Output:** The written path, sample count and column names.
- **Example code:**

```python
from dnakit.sequence.skbio_p20_io_functions import read_skbio_file, write_skbio_sample_metadata
out = write_skbio_sample_metadata({"group": ["x", "y"], "depth": [10, 20]}, ["a", "b"], "meta.tsv", verbose=False)
print(out["columns"])
print(open("meta.tsv").read())
back = read_skbio_file("meta.tsv", "sample_metadata", verbose=False)
print(back["column_types"])
```

- **Example output:**

```text
['group', 'depth']
id	group	depth
#sk:types	categorical	numeric
a	x	10
b	y	20

{'group': 'categorical', 'depth': 'numeric'}
```

## 5) `dnakit.sequence.skbio_p20_io_functions.write_skbio_biom_table()`

- **Purpose:** Write a count table as a BIOM (HDF5) file.
- **Method:** Builds a `biom.Table` (features by samples) from the sample-by-feature counts and writes it with scikit-bio's `biom` writer. Needs the `biom-format` package.
- **API:** `dnakit.sequence.skbio_p20_io_functions.write_skbio_biom_table(counts, output_path, *, sample_ids=None, feature_ids=None)`.
- **Input:** Counts (one row per sample), optional sample and feature ids, an output path.
- **Output:** The written path, sample count and feature count.
- **Example code:**

```python
from dnakit.sequence.skbio_p20_io_functions import read_skbio_file, write_skbio_biom_table
out = write_skbio_biom_table([[1, 2], [3, 4]], "table.biom", sample_ids=["s1", "s2"], feature_ids=["f1", "f2"], verbose=False)
print(out["sample_count"], out["feature_count"])
back = read_skbio_file("table.biom", "biom", verbose=False)
print(back["sample_ids"], back["feature_ids"], back["counts_features_by_samples"])
```

- **Example output:**

```text
2 2
['s1', 's2'] ['f1', 'f2'] [[1.0, 3.0], [2.0, 4.0]]
```

## 6) `dnakit.sequence.skbio_p20_io_functions.lowercase_masked_sequence()`

- **Purpose:** Render a sequence with chosen intervals in lower case (soft-masking).
- **Method:** Builds a mask from the `[start, end)` intervals and calls `skbio.Sequence.lowercase`.
- **API:** `dnakit.sequence.skbio_p20_io_functions.lowercase_masked_sequence(sequence, mask_intervals, *, sequence_type='DNA')`.
- **Input:** An upper-case IUPAC sequence and a list of 0-based half-open intervals.
- **Output:** The masked `sequence`, the number of masked positions and the length.
- **Example code:**

```python
from dnakit.sequence.skbio_p20_io_functions import lowercase_masked_sequence
r = lowercase_masked_sequence("ACGTACGT", [[1, 3], [6, 8]], verbose=False)
print(r["sequence"], r["masked_positions"])
```

- **Example output:**

```text
AcgTACgt 4
```

## 7) `dnakit.sequence.skbio_p20_io_functions.find_subsequence()`

- **Purpose:** Find the first occurrence of a subsequence.
- **Method:** Calls `skbio.Sequence.index`; `start` and `end` bound the search. An absent subsequence gives `found=False` instead of an error.
- **API:** `dnakit.sequence.skbio_p20_io_functions.find_subsequence(sequence, subsequence, *, start=None, end=None, sequence_type='DNA')`.
- **Input:** A sequence, a subsequence, optional bounds.
- **Output:** `found` and the 0-based `index`.
- **Example code:**

```python
from dnakit.sequence.skbio_p20_io_functions import find_subsequence
print(find_subsequence("ACGTACGT", "GT", verbose=False)["index"])
print(find_subsequence("ACGTACGT", "GT", start=4, verbose=False)["index"])
print(find_subsequence("ACGT", "TT", verbose=False)["found"])
```

- **Example output:**

```text
2
6
False
```

## 8) `dnakit.sequence.skbio_p20_io_functions.split_contiguous_segments()`

- **Purpose:** Split a sequence into contiguous segments selected by a mask.
- **Method:** Wraps `skbio.Sequence.iter_contiguous`. The mask is a list of intervals or every position whose character is in `included_characters`; `invert` selects the complement and `min_length` drops short segments.
- **API:** `dnakit.sequence.skbio_p20_io_functions.split_contiguous_segments(sequence, *, included_intervals=None, included_characters=None, min_length=1, invert=False, sequence_type='DNA')`.
- **Input:** A sequence and exactly one of `included_intervals` and `included_characters`.
- **Output:** `segments` with `start`, `end` and `sequence`.
- **Example code:**

```python
from dnakit.sequence.skbio_p20_io_functions import split_contiguous_segments
r = split_contiguous_segments("ACG-TNA", included_characters="-N", invert=True, verbose=False)
for s in r["segments"]:
    print(s["start"], s["end"], s["sequence"])
```

- **Example output:**

```text
0 3 ACG
4 5 T
6 7 A
```

## 9) `dnakit.sequence.skbio_p20_io_functions.locate_sequence_positions()`

- **Purpose:** Locate gap, degenerate or definite positions of an IUPAC sequence.
- **Method:** `kind` selects `GrammaredSequence.gaps`, `degenerates` or `definites`; the result also states whether the sequence has any definite character (`has_definites`). One function covers four rows of the work list.
- **API:** `dnakit.sequence.skbio_p20_io_functions.locate_sequence_positions(sequence, kind, *, sequence_type='DNA')`.
- **Input:** A sequence and `kind`.
- **Output:** `positions`, `runs` as `[start, end)`, `count`, `length`, `has_definites`.
- **Example code:**

```python
from dnakit.sequence.skbio_p20_io_functions import locate_sequence_positions
for kind in ("gaps", "degenerates", "definites"):
    r = locate_sequence_positions("ACG-TNA", kind, verbose=False)
    print(kind, r["positions"], r["runs"], r["has_definites"])
print(locate_sequence_positions("NNN", "definites", verbose=False)["has_definites"])
```

- **Example output:**

```text
gaps [3] [[3, 4]] True
degenerates [5] [[5, 6]] True
definites [0, 1, 2, 4, 6] [[0, 3], [4, 5], [6, 7]] True
False
```

## 10) `dnakit.sequence.skbio_p20_io_functions.find_nucleotide_motifs()`

- **Purpose:** Find purine or pyrimidine runs in a nucleotide sequence.
- **Method:** Calls `GrammaredSequence.find_motifs` with `purine-run` or `pyrimidine-run`; `ignore_intervals` skips positions and `min_length` filters short runs.
- **API:** `dnakit.sequence.skbio_p20_io_functions.find_nucleotide_motifs(sequence, motif_type, *, min_length=1, ignore_intervals=None, sequence_type='DNA')`.
- **Input:** A sequence, the motif type and optional filters.
- **Output:** `motifs` with `start`, `end` and `sequence`, and `count`.
- **Example code:**

```python
from dnakit.sequence.skbio_p20_io_functions import find_nucleotide_motifs
r = find_nucleotide_motifs("ACGGGGAGGCGGAG", "purine-run", min_length=2, verbose=False)
for m in r["motifs"]:
    print(m["start"], m["end"], m["sequence"])
print(find_nucleotide_motifs("ACGGGGAGGCGGAG", "pyrimidine-run", min_length=2, verbose=False)["count"])
```

- **Example output:**

```text
2 9 GGGGAGG
10 14 GGAG
0
```

## 11) `dnakit.sequence.skbio_p20_io_functions.is_reverse_complement_pair()`

- **Purpose:** Test whether two sequences are reverse complements of each other.
- **Method:** Calls `NucleotideMixin.is_reverse_complement`, which handles IUPAC degenerate codes.
- **API:** `dnakit.sequence.skbio_p20_io_functions.is_reverse_complement_pair(sequence_a, sequence_b, *, sequence_type='DNA')`.
- **Input:** Two sequences of the same type.
- **Output:** `is_reverse_complement` (true or false).
- **Example code:**

```python
from dnakit.sequence.skbio_p20_io_functions import is_reverse_complement_pair
print(is_reverse_complement_pair("AACG", "CGTT", verbose=False)["is_reverse_complement"])
print(is_reverse_complement_pair("AACG", "ACGT", verbose=False)["is_reverse_complement"])
```

- **Example output:**

```text
True
False
```

## 12) `dnakit.sequence.skbio_p20_io_functions.msa_consensus_sequence()`

- **Purpose:** Majority consensus of a multiple sequence alignment.
- **Method:** Builds a `TabularMSA` and calls `consensus()`; ties follow scikit-bio's rule.
- **API:** `dnakit.sequence.skbio_p20_io_functions.msa_consensus_sequence(alignment, *, sequence_type='DNA')`.
- **Input:** A mapping of name to aligned sequence (equal lengths).
- **Output:** `consensus`, `sequence_count`, `alignment_length`.
- **Example code:**

```python
from dnakit.sequence.skbio_p20_io_functions import msa_consensus_sequence
aln = {"a": "AC-T", "b": "ACGT", "c": "ACGA"}
r = msa_consensus_sequence(aln, verbose=False)
print(r["consensus"], r["sequence_count"], r["alignment_length"])
```

- **Example output:**

```text
ACGT 3 4
```

## 13) `dnakit.sequence.skbio_p20_io_functions.msa_conservation()`

- **Purpose:** Per-column conservation of an alignment.
- **Method:** Calls `TabularMSA.conservation` with scikit-bio defaults (inverse Shannon uncertainty, columns with gaps are null). `gap_mode` and `degenerate_mode` change the handling.
- **API:** `dnakit.sequence.skbio_p20_io_functions.msa_conservation(alignment, *, metric='inverse_shannon_uncertainty', degenerate_mode='error', gap_mode='nan', sequence_type='DNA')`.
- **Input:** An aligned mapping and the options above.
- **Output:** `conservation`, one value per column.
- **Example code:**

```python
from dnakit.sequence.skbio_p20_io_functions import msa_conservation
aln = {"a": "AC-T", "b": "ACGT", "c": "ACGA"}
print(msa_conservation(aln, verbose=False)["conservation"])
print(msa_conservation(aln, gap_mode="ignore", verbose=False)["conservation"])
```

- **Example output:**

```text
[1.0, 1.0, None, 0.5408520829727552]
[1.0, 1.0, 1.0, 0.5408520829727552]
```

## 14) `dnakit.sequence.skbio_p20_io_functions.msa_gap_frequencies()`

- **Purpose:** Gap counts (or fractions) per alignment column or per sequence.
- **Method:** Calls `TabularMSA.gap_frequencies`; `axis='sequence'` gives one value per column, `axis='position'` one per sequence; `relative` returns fractions.
- **API:** `dnakit.sequence.skbio_p20_io_functions.msa_gap_frequencies(alignment, *, axis='sequence', relative=False, sequence_type='DNA')`.
- **Input:** An aligned mapping, `axis`, `relative`.
- **Output:** `frequencies` (and `names` for the per-sequence axis).
- **Example code:**

```python
from dnakit.sequence.skbio_p20_io_functions import msa_gap_frequencies
aln = {"a": "AC-T", "b": "ACGT", "c": "ACGA"}
print(msa_gap_frequencies(aln, verbose=False)["frequencies"])
r = msa_gap_frequencies(aln, axis="position", relative=True, verbose=False)
print(r["names"], r["frequencies"])
```

- **Example output:**

```text
[0, 0, 1, 0]
['a', 'b', 'c'] [0.25, 0.0, 0.0]
```

## 15) `dnakit.sequence.skbio_p20_io_functions.parse_cigar_alignment_path()`

- **Purpose:** Parse a CIGAR string into a pairwise alignment path.
- **Method:** Calls `PairAlignPath.from_cigar` and reports the run lengths, states, start offsets, alignment length and the CIGAR re-encoded from the path.
- **API:** `dnakit.sequence.skbio_p20_io_functions.parse_cigar_alignment_path(cigar, *, starts=None)`.
- **Input:** A CIGAR string and optional `starts`.
- **Output:** `lengths`, `states`, `starts`, `alignment_length`, `cigar`.
- **Example code:**

```python
from dnakit.sequence.skbio_p20_io_functions import parse_cigar_alignment_path
r = parse_cigar_alignment_path("3M1D2M", verbose=False)
print(r["lengths"], r["states"], r["alignment_length"], r["cigar"])
```

- **Example output:**

```text
[3, 1, 2] [[0, 2, 0]] 6 3M1D2M
```

## 16) `dnakit.sequence.skbio_p20_io_functions.apply_alignment_path()`

- **Purpose:** Rebuild the gapped alignment of the original sequences from an alignment path.
- **Method:** Builds an `AlignPath` from a pairwise `cigar`, or from run `lengths` and per-run 0/1 gap flags (`states`) for any number of sequences, and calls `AlignPath.to_aligned` with `gap_char`.
- **API:** `dnakit.sequence.skbio_p20_io_functions.apply_alignment_path(sequences, *, cigar=None, lengths=None, states=None, starts=None, gap_char='-', sequence_type='DNA')`.
- **Input:** A mapping of name to unaligned sequence and either `cigar` or `lengths` with `states` and `starts`.
- **Output:** `alignment` (name to gapped sequence) and `alignment_length`.
- **Example code:**

```python
from dnakit.sequence.skbio_p20_io_functions import apply_alignment_path
seqs = {"x": "ACGTAC", "y": "ACGAC"}
print(apply_alignment_path(seqs, cigar="3M1D2M", verbose=False)["alignment"])
r = apply_alignment_path(seqs, lengths=[3, 1, 2], states=[[0, 0], [0, 1], [0, 0]], starts=[0, 0], verbose=False)
print(r["alignment"])
```

- **Example output:**

```text
{'x': 'ACGTAC', 'y': 'ACG-AC'}
{'x': 'ACGTAC', 'y': 'ACG-AC'}
```
