# AUGUSTUS: hints, alignment filters and coverage tracks

The functions in `dnakit.annotation` listed here are typed facades over AUGUSTUS 3.5.0 programs and helper scripts ([AUGUSTUS](https://github.com/Gaius-Augustus/Augustus)). Each one runs a single fixed program with an explicit argument list (no shell). Install AUGUSTUS (`conda install -c bioconda augustus`) and put its `bin` folder on PATH, or pass `augustus_bin_dir`. If the program is missing the function raises `BackendUnavailableError` with the install hint. The species parameter sets are read from `AUGUSTUS_CONFIG_PATH` (or `config_path`); training functions write there, so the folder must be writable. Output paths must be new; nothing is overwritten. Each result is a dictionary with `status`, `backend`, `version`, `command`, `return_code`, `outputs`, `log_tail` and `provenance`, plus small parsed counts.

## 1) `dnakit.annotation.augustus_alignments_to_wiggle()`

- **Purpose:** Convert PSL or SHRIMP alignments to wiggle coverage (``aln2wig -f``); ``span`` uses span notation.
- **Method:** Runs `span` through the DNAKit bounded executor.
- **API:** `dnakit.annotation.augustus_alignments_to_wiggle(alignment_path, output_path, *, span=False, track_name=None, augustus_bin_dir=None, verbose=True)`.
- **Input:** Required: `alignment_path`, `output_path`. Optional keyword-only: `span`, `track_name`, `augustus_bin_dir` (including `augustus_bin_dir`, `config_path` where applicable). Input files must exist; output paths must not.
- **Output:** A dictionary with `status`, `backend`, `version`, `command`, `return_code`, `outputs`, `log_tail`, `provenance` and parsed counts such as `lines`, `records`, `genes` or `transcripts`.
- **Example code:**

```python
import random, tempfile
from pathlib import Path
from dnakit.annotation import augustus_p04_hints_functions as m
w = Path(tempfile.mkdtemp())
(w / "a.psl").write_text("100\t0\t0\t0\t0\t0\t1\t300\t+\tq1\t100\t0\t100\tchr1\t4000\t1000\t1400\t2\t50,50,\t0,50,\t1000,1350,\n")
r = m.augustus_alignments_to_wiggle(str(w / "a.psl"), str(w / "aln.wig"), verbose=False)
print(r["lines"], r["preview"][1])
```

- **Example output:**

```text
102 variableStep chrom=chr1
```

## 2) `dnakit.annotation.augustus_bam_to_hints()`

- **Purpose:** Derive intron (and optionally exon) hints from spliced RNA-Seq alignments in BAM (``bam2hints``).
- **Method:** Runs `bam2hints`, `exon_hints`, `no_multiplicity` through the DNAKit bounded executor. The BAM must be sorted by target name and start. ``exon_hints`` adds exonpart, exon and splice-site hints; ``no_multiplicity`` keeps identical intron hints separate instead of summarising them with ``mult=n``. Unset options keep the program default (priority 4, max gap 14, min intron 32, max intron 350000, source E).
- **API:** `dnakit.annotation.augustus_bam_to_hints(bam_path, output_path, *, priority=None, max_gap_length=None, min_intron_length=None, max_intron_length=None, min_end_block_length=None, max_query_gap_length=None, exon_hints=False, ep_cutoff=None, source=None, no_multiplicity=False, remove_redundant=False, max_coverage=None, splice_site_hints=False, truncated_splice_site_hints=False, score=None, max_gene_length=None, augustus_bin_dir=None, verbose=True)`.
- **Input:** Required: `bam_path`, `output_path`. Optional keyword-only: `priority`, `max_gap_length`, `min_intron_length`, `max_intron_length`, `min_end_block_length`, `max_query_gap_length`, `exon_hints`, `ep_cutoff`, `source`, `no_multiplicity`, `remove_redundant`, `max_coverage`, `splice_site_hints`, `truncated_splice_site_hints`, `score`, `max_gene_length`, `augustus_bin_dir` (including `augustus_bin_dir`, `config_path` where applicable). Input files must exist; output paths must not.
- **Output:** A dictionary with `status`, `backend`, `version`, `command`, `return_code`, `outputs`, `log_tail`, `provenance` and parsed counts such as `lines`, `records`, `genes` or `transcripts`.
- **Example code:** not shown; the example needs data that is too large or an external tool that is not part of the AUGUSTUS package (see the tests and the manifest notes).

## 3) `dnakit.annotation.augustus_bam_to_wiggle()`

- **Purpose:** Convert BAM alignments to wiggle coverage for exonpart hints.
- **Method:** Runs `faToTwoBit`, `genome_path`, `region`, `track_name` through the DNAKit bounded executor. ``method="bam2wig"`` runs the AUGUSTUS program (sorted BAM; ``region`` such as ``chr1:10-250`` needs a BAM index; ``track_name``). ``method="bamToWig.py"`` runs the Python alternative, which also needs ``genome_path`` (FASTA) and samtools plus the UCSC tools ``twoBitInfo`` and ``faToTwoBit`` on PATH.
- **API:** `dnakit.annotation.augustus_bam_to_wiggle(bam_path, output_path, *, method='bam2wig', genome_path=None, track_name=None, region=None, augustus_bin_dir=None, verbose=True)`.
- **Input:** Required: `bam_path`, `output_path`. Optional keyword-only: `method`, `genome_path`, `track_name`, `region`, `augustus_bin_dir` (including `augustus_bin_dir`, `config_path` where applicable). Input files must exist; output paths must not.
- **Output:** A dictionary with `status`, `backend`, `version`, `command`, `return_code`, `outputs`, `log_tail`, `provenance` and parsed counts such as `lines`, `records`, `genes` or `transcripts`.
- **Example code:** not shown; the example needs data that is too large or an external tool that is not part of the AUGUSTUS package (see the tests and the manifest notes).

## 4) `dnakit.annotation.augustus_bedgraph_to_wiggle()`

- **Purpose:** Convert bedGraph coverage to variableStep wiggle (``bedgraph2wig.pl``); ``span`` uses span notation.
- **Method:** Runs `bedgraph2wig.pl`, `span` through the DNAKit bounded executor.
- **API:** `dnakit.annotation.augustus_bedgraph_to_wiggle(bedgraph_path, output_path, *, span=False, augustus_bin_dir=None, verbose=True)`.
- **Input:** Required: `bedgraph_path`, `output_path`. Optional keyword-only: `span`, `augustus_bin_dir` (including `augustus_bin_dir`, `config_path` where applicable). Input files must exist; output paths must not.
- **Output:** A dictionary with `status`, `backend`, `version`, `command`, `return_code`, `outputs`, `log_tail`, `provenance` and parsed counts such as `lines`, `records`, `genes` or `transcripts`.
- **Example code:**

```python
import random, tempfile
from pathlib import Path
from dnakit.annotation import augustus_p04_hints_functions as m
w = Path(tempfile.mkdtemp())
(w / "c.bedgraph").write_text("chr1\t0\t100\t3\nchr1\t100\t200\t5\n")
r = m.augustus_bedgraph_to_wiggle(str(w / "c.bedgraph"), str(w / "c.wig"), verbose=False)
print(r["preview"])
```

- **Example output:**

```text
['track type=wiggle_0 name=chr1\t0\t100\t3 ', 'chr1\t100\t200\t5']
```

## 5) `dnakit.annotation.augustus_choose_wiggle_intervals()`

- **Purpose:** Choose, per target sequence, an interval with good wiggle coverage and write it as BED (``wigchoose.pl``).
- **Method:** Runs `wigchoose.pl` through the DNAKit bounded executor.
- **API:** `dnakit.annotation.augustus_choose_wiggle_intervals(wiggle_path, output_path, *, min_coverage=None, min_relative_coverage=None, max_gap=None, augustus_bin_dir=None, verbose=True)`.
- **Input:** Required: `wiggle_path`, `output_path`. Optional keyword-only: `min_coverage`, `min_relative_coverage`, `max_gap`, `augustus_bin_dir` (including `augustus_bin_dir`, `config_path` where applicable). Input files must exist; output paths must not.
- **Output:** A dictionary with `status`, `backend`, `version`, `command`, `return_code`, `outputs`, `log_tail`, `provenance` and parsed counts such as `lines`, `records`, `genes` or `transcripts`.
- **Example code:**

```python
import random, tempfile
from pathlib import Path
from dnakit.annotation import augustus_p04_hints_functions as m
w = Path(tempfile.mkdtemp())
(w / "cov.wig").write_text("variableStep chrom=chr1\n" + "".join(f"{p} 5\n" for p in range(101, 260)))
r = m.augustus_choose_wiggle_intervals(str(w / "cov.wig"), str(w / "choose.bed"), verbose=False)
print(r["preview"])
```

- **Example output:**

```text
['chr1\t101\t259\t795']
```

## 6) `dnakit.annotation.augustus_compile_splice_candidates()`

- **Purpose:** Compile intron candidates from exon-candidate and coverage hints (``compileSpliceCands``); GFF goes to ``output_path``.
- **Method:** Runs `compileSpliceCands`, `output_path` through the DNAKit bounded executor.
- **API:** `dnakit.annotation.augustus_compile_splice_candidates(candidates_path, genome_path, output_path, *, compare_count=None, coverage_tolerance=None, margin=None, max_intron_length=None, augustus_bin_dir=None, verbose=True)`.
- **Input:** Required: `candidates_path`, `genome_path`, `output_path`. Optional keyword-only: `compare_count`, `coverage_tolerance`, `margin`, `max_intron_length`, `augustus_bin_dir` (including `augustus_bin_dir`, `config_path` where applicable). Input files must exist; output paths must not.
- **Output:** A dictionary with `status`, `backend`, `version`, `command`, `return_code`, `outputs`, `log_tail`, `provenance` and parsed counts such as `lines`, `records`, `genes` or `transcripts`.
- **Example code:** not shown; the example needs data that is too large or an external tool that is not part of the AUGUSTUS package (see the tests and the manifest notes).

## 7) `dnakit.annotation.augustus_exonerate_to_hints()`

- **Purpose:** Convert Exonerate protein2genome GFF alignments to AUGUSTUS hints (``exonerate2hints.pl``).
- **Method:** Runs `exonerate2hints.pl` through the DNAKit bounded executor.
- **API:** `dnakit.annotation.augustus_exonerate_to_hints(exonerate_path, output_path, *, priority=None, min_intron_length=None, max_intron_length=None, cds_part_cutoff=None, source=None, splice_site_hints=False, augustus_bin_dir=None, verbose=True)`.
- **Input:** Required: `exonerate_path`, `output_path`. Optional keyword-only: `priority`, `min_intron_length`, `max_intron_length`, `cds_part_cutoff`, `source`, `splice_site_hints`, `augustus_bin_dir` (including `augustus_bin_dir`, `config_path` where applicable). Input files must exist; output paths must not.
- **Output:** A dictionary with `status`, `backend`, `version`, `command`, `return_code`, `outputs`, `log_tail`, `provenance` and parsed counts such as `lines`, `records`, `genes` or `transcripts`.
- **Example code:** not shown; the example needs data that is too large or an external tool that is not part of the AUGUSTUS package (see the tests and the manifest notes).

## 8) `dnakit.annotation.augustus_exoniphy_to_hints()`

- **Purpose:** Convert exoniphy exon predictions (tab table: name, chrom, strand, cdsStart, cdsEnd) to hints (``exoniphyDb2hints.pl``).
- **Method:** Runs `exoniphyDb2hints.pl` through the DNAKit bounded executor.
- **API:** `dnakit.annotation.augustus_exoniphy_to_hints(table_path, output_path, *, priority=None, cutoff=None, augustus_bin_dir=None, verbose=True)`.
- **Input:** Required: `table_path`, `output_path`. Optional keyword-only: `priority`, `cutoff`, `augustus_bin_dir` (including `augustus_bin_dir`, `config_path` where applicable). Input files must exist; output paths must not.
- **Output:** A dictionary with `status`, `backend`, `version`, `command`, `return_code`, `outputs`, `log_tail`, `provenance` and parsed counts such as `lines`, `records`, `genes` or `transcripts`.
- **Example code:** not shown; the example needs data that is too large or an external tool that is not part of the AUGUSTUS package (see the tests and the manifest notes).

## 9) `dnakit.annotation.augustus_filter_bam()`

- **Purpose:** Filter BAM alignments by identity, coverage, uniqueness and pairing (``filterBam``).
- **Method:** Runs `common_gene_path`, `filterBam`, `pair_bed_path`, `samtools` through the DNAKit bounded executor. The BAM must be sorted by query name. Program defaults: min identity 92, min coverage 80, unique threshold 0.96. ``pair_bed_path`` and ``common_gene_path`` are new optional output files. Returns the record count of the filtered BAM when ``samtools`` is available.
- **API:** `dnakit.annotation.augustus_filter_bam(bam_path, output_path, *, unique=False, unique_threshold=None, best=False, min_identity=None, min_coverage=None, no_introns=False, insert_limit=None, paired=False, max_intron_length=None, pair_bed_path=None, common_gene_path=None, threads=None, augustus_bin_dir=None, verbose=True)`.
- **Input:** Required: `bam_path`, `output_path`. Optional keyword-only: `unique`, `unique_threshold`, `best`, `min_identity`, `min_coverage`, `no_introns`, `insert_limit`, `paired`, `max_intron_length`, `pair_bed_path`, `common_gene_path`, `threads`, `augustus_bin_dir` (including `augustus_bin_dir`, `config_path` where applicable). Input files must exist; output paths must not.
- **Output:** A dictionary with `status`, `backend`, `version`, `command`, `return_code`, `outputs`, `log_tail`, `provenance` and parsed counts such as `lines`, `records`, `genes` or `transcripts`.
- **Example code:** not shown; the example needs data that is too large or an external tool that is not part of the AUGUSTUS package (see the tests and the manifest notes).

## 10) `dnakit.annotation.augustus_filter_maf()`

- **Purpose:** Keep MAF blocks that intersect ``interval`` (``seqID:start-end``) or restrict them to ``species`` (``filterMaf.pl``).
- **Method:** Runs `filterMaf.pl`, `interval`, `species` through the DNAKit bounded executor.
- **API:** `dnakit.annotation.augustus_filter_maf(maf_path, output_path, *, species=None, min_sequences=None, interval=None, augustus_bin_dir=None, verbose=True)`.
- **Input:** Required: `maf_path`, `output_path`. Optional keyword-only: `species`, `min_sequences`, `interval`, `augustus_bin_dir` (including `augustus_bin_dir`, `config_path` where applicable). Input files must exist; output paths must not.
- **Output:** A dictionary with `status`, `backend`, `version`, `command`, `return_code`, `outputs`, `log_tail`, `provenance` and parsed counts such as `lines`, `records`, `genes` or `transcripts`.
- **Example code:**

```python
import random, tempfile
from pathlib import Path
from dnakit.annotation import augustus_p04_hints_functions as m
w = Path(tempfile.mkdtemp())
(w / "t.maf").write_text("##maf version=1\n\na score=100\ns hg.chr1 100 20 + 1000 ACGT-ACGTACGTACGTACGT\ns mm.chr2 200 20 + 1000 ACGTTACGTACGTACGTACGT\n\n")
inside = m.augustus_filter_maf(str(w / "t.maf"), str(w / "in.maf"), interval="hg.chr1:50-300", verbose=False)
outside = m.augustus_filter_maf(str(w / "t.maf"), str(w / "out.maf"), interval="hg.chr1:500-600", verbose=False)
print(inside["lines"], outside["lines"])
```

- **Example output:**

```text
3 0
```

## 11) `dnakit.annotation.augustus_filter_psl()`

- **Purpose:** Filter a PSL alignment file (BLAT, GMAP) by identity and coverage (``filterPSL.pl``); the PSL must be sorted by query name.
- **Method:** Runs `filterPSL.pl` through the DNAKit bounded executor.
- **API:** `dnakit.annotation.augustus_filter_psl(psl_path, output_path, *, min_identity=None, min_coverage=None, unique=False, unique_threshold=None, best=False, no_introns=False, paired=False, max_intron_length=None, pair_bed_path=None, common_gene_path=None, augustus_bin_dir=None, verbose=True)`.
- **Input:** Required: `psl_path`, `output_path`. Optional keyword-only: `min_identity`, `min_coverage`, `unique`, `unique_threshold`, `best`, `no_introns`, `paired`, `max_intron_length`, `pair_bed_path`, `common_gene_path`, `augustus_bin_dir` (including `augustus_bin_dir`, `config_path` where applicable). Input files must exist; output paths must not.
- **Output:** A dictionary with `status`, `backend`, `version`, `command`, `return_code`, `outputs`, `log_tail`, `provenance` and parsed counts such as `lines`, `records`, `genes` or `transcripts`.
- **Example code:**

```python
import random, tempfile
from pathlib import Path
from dnakit.annotation import augustus_p04_hints_functions as m
w = Path(tempfile.mkdtemp())
(w / "a.psl").write_text("100\t0\t0\t0\t0\t0\t1\t300\t+\tq1\t100\t0\t100\tchr1\t4000\t1000\t1400\t2\t50,50,\t0,50,\t1000,1350,\n")
keep = m.augustus_filter_psl(str(w / "a.psl"), str(w / "keep.psl"), verbose=False)
drop = m.augustus_filter_psl(str(w / "a.psl"), str(w / "drop.psl"), no_introns=True, verbose=False)
print(keep["lines"], drop["lines"])
```

- **Example output:**

```text
1 0
```

## 12) `dnakit.annotation.augustus_filter_shrimp()`

- **Purpose:** Filter a SHRIMP alignment file by score (``filterShrimp.pl``); the file must be sorted by query name.
- **Method:** Runs `filterShrimp.pl` through the DNAKit bounded executor.
- **API:** `dnakit.annotation.augustus_filter_shrimp(shrimp_path, output_path, *, min_score=None, unique=False, unique_threshold=None, best=False, common_gene_path=None, augustus_bin_dir=None, verbose=True)`.
- **Input:** Required: `shrimp_path`, `output_path`. Optional keyword-only: `min_score`, `unique`, `unique_threshold`, `best`, `common_gene_path`, `augustus_bin_dir` (including `augustus_bin_dir`, `config_path` where applicable). Input files must exist; output paths must not.
- **Output:** A dictionary with `status`, `backend`, `version`, `command`, `return_code`, `outputs`, `log_tail`, `provenance` and parsed counts such as `lines`, `records`, `genes` or `transcripts`.
- **Example code:** not shown; the example needs data that is too large or an external tool that is not part of the AUGUSTUS package (see the tests and the manifest notes).

## 13) `dnakit.annotation.augustus_filter_splice_hints()`

- **Purpose:** Keep intron hints whose splice sites match ``splice_pattern`` (usually GTAG, sometimes GCAG) in the genome (``filterSpliceHints.pl``).
- **Method:** Runs `filterSpliceHints.pl`, `splice_pattern` through the DNAKit bounded executor.
- **API:** `dnakit.annotation.augustus_filter_splice_hints(genome_path, hints_path, output_path, *, splice_pattern='GTAG', augustus_bin_dir=None, verbose=True)`.
- **Input:** Required: `genome_path`, `hints_path`, `output_path`. Optional keyword-only: `splice_pattern`, `augustus_bin_dir` (including `augustus_bin_dir`, `config_path` where applicable). Input files must exist; output paths must not.
- **Output:** A dictionary with `status`, `backend`, `version`, `command`, `return_code`, `outputs`, `log_tail`, `provenance` and parsed counts such as `lines`, `records`, `genes` or `transcripts`.
- **Example code:**

```python
import random, tempfile
from pathlib import Path
from dnakit.annotation import augustus_p04_hints_functions as m
w = Path(tempfile.mkdtemp())
rng = random.Random(7)
(w / "genome.fa").write_text(">chr1\n" + "".join(rng.choices("ACGT", k=4000)) + "\n")
(w / "h.gff").write_text("chr1\tb2h\tintron\t1051\t1350\t0\t+\t.\tgrp=a;pri=4;src=E\nchr1\tb2h\tintron\t1051\t1350\t0\t+\t.\tgrp=b;pri=4;src=E\n")
r = m.augustus_filter_splice_hints(str(w / "genome.fa"), str(w / "h.gff"), str(w / "splice.gff"), verbose=False)
print(r["status"], r["lines"])
```

- **Example output:**

```text
completed 0
```

## 14) `dnakit.annotation.augustus_join_multiplicity_hints()`

- **Purpose:** Collapse identical hints of one source into one hint with ``mult=n`` (``join_mult_hints.pl``); the hints must be sorted.
- **Method:** Runs `join_mult_hints.pl` through the DNAKit bounded executor.
- **API:** `dnakit.annotation.augustus_join_multiplicity_hints(hints_path, output_path, *, augustus_bin_dir=None, verbose=True)`.
- **Input:** Required: `hints_path`, `output_path`. Optional keyword-only: `augustus_bin_dir` (including `augustus_bin_dir`, `config_path` where applicable). Input files must exist; output paths must not.
- **Output:** A dictionary with `status`, `backend`, `version`, `command`, `return_code`, `outputs`, `log_tail`, `provenance` and parsed counts such as `lines`, `records`, `genes` or `transcripts`.
- **Example code:**

```python
import random, tempfile
from pathlib import Path
from dnakit.annotation import augustus_p04_hints_functions as m
w = Path(tempfile.mkdtemp())
(w / "h.gff").write_text("chr1\tb2h\tintron\t1051\t1350\t0\t+\t.\tgrp=a;pri=4;src=E\nchr1\tb2h\tintron\t1051\t1350\t0\t+\t.\tgrp=b;pri=4;src=E\n")
r = m.augustus_join_multiplicity_hints(str(w / "h.gff"), str(w / "joined.gff"), verbose=False)
print(r["lines"], "mult=2" in r["preview"][0])
```

- **Example output:**

```text
1 True
```

## 15) `dnakit.annotation.augustus_make_utr_training_set()`

- **Purpose:** Build a UTR training set (GenBank and GFF) from start/stop codons and filtered EST PSL alignments (``makeUtrTrainingSet.pl``).
- **Method:** Runs `makeUtrTrainingSet.pl` through the DNAKit bounded executor. Writes ``<output_prefix>.gb`` and ``<output_prefix>.gff`` (both must not exist).
- **API:** `dnakit.annotation.augustus_make_utr_training_set(codons_gff_path, genome_path, est_psl_path, output_prefix, *, only_best=False, distance=None, augustus_bin_dir=None, verbose=True)`.
- **Input:** Required: `codons_gff_path`, `genome_path`, `est_psl_path`, `output_prefix`. Optional keyword-only: `only_best`, `distance`, `augustus_bin_dir` (including `augustus_bin_dir`, `config_path` where applicable). Input files must exist; output paths must not.
- **Output:** A dictionary with `status`, `backend`, `version`, `command`, `return_code`, `outputs`, `log_tail`, `provenance` and parsed counts such as `lines`, `records`, `genes` or `transcripts`.
- **Example code:** not shown; the example needs data that is too large or an external tool that is not part of the AUGUSTUS package (see the tests and the manifest notes).

## 16) `dnakit.annotation.augustus_map_junction_sam()`

- **Purpose:** Map a SAM of reads aligned to an exon-exon junction file back to genome coordinates (``samMap.pl``); the default flank is 75 bp.
- **Method:** Runs `samMap.pl` through the DNAKit bounded executor.
- **API:** `dnakit.annotation.augustus_map_junction_sam(sam_path, map_psl_path, output_path, *, flank=None, augustus_bin_dir=None, verbose=True)`.
- **Input:** Required: `sam_path`, `map_psl_path`, `output_path`. Optional keyword-only: `flank`, `augustus_bin_dir` (including `augustus_bin_dir`, `config_path` where applicable). Input files must exist; output paths must not.
- **Output:** A dictionary with `status`, `backend`, `version`, `command`, `return_code`, `outputs`, `log_tail`, `provenance` and parsed counts such as `lines`, `records`, `genes` or `transcripts`.
- **Example code:** not shown; the example needs data that is too large or an external tool that is not part of the AUGUSTUS package (see the tests and the manifest notes).

## 17) `dnakit.annotation.augustus_maf_to_conservation_wiggle()`

- **Purpose:** Convert a MAF alignment to a conservation wiggle track (``maf2conswig.pl``).
- **Method:** Runs `maf2conswig.pl` through the DNAKit bounded executor.
- **API:** `dnakit.annotation.augustus_maf_to_conservation_wiggle(maf_path, output_path, *, window_size=None, augustus_bin_dir=None, verbose=True)`.
- **Input:** Required: `maf_path`, `output_path`. Optional keyword-only: `window_size`, `augustus_bin_dir` (including `augustus_bin_dir`, `config_path` where applicable). Input files must exist; output paths must not.
- **Output:** A dictionary with `status`, `backend`, `version`, `command`, `return_code`, `outputs`, `log_tail`, `provenance` and parsed counts such as `lines`, `records`, `genes` or `transcripts`.
- **Example code:**

```python
import random, tempfile
from pathlib import Path
from dnakit.annotation import augustus_p04_hints_functions as m
w = Path(tempfile.mkdtemp())
(w / "t.maf").write_text("##maf version=1\n\na score=100\ns hg.chr1 100 20 + 1000 ACGT-ACGTACGTACGTACGT\ns mm.chr2 200 20 + 1000 ACGTTACGTACGTACGTACGT\n\n")
r = m.augustus_maf_to_conservation_wiggle(str(w / "t.maf"), str(w / "cons.wig"), verbose=False)
print(r["preview"][0])
```

- **Example output:**

```text
fixedStep chrom=chr1 start=101 step=1
```

## 18) `dnakit.annotation.augustus_peptides_to_hints()`

- **Purpose:** Convert peptide-to-genome PSL alignments (with peptide GFF and FASTA) to hints (``peptides2hints.pl``).
- **Method:** Runs `peptides2hints.pl` through the DNAKit bounded executor.
- **API:** `dnakit.annotation.augustus_peptides_to_hints(psl_path, gff_path, fasta_path, output_path, *, source='P', priority=4, augustus_bin_dir=None, verbose=True)`.
- **Input:** Required: `psl_path`, `gff_path`, `fasta_path`, `output_path`. Optional keyword-only: `source`, `priority`, `augustus_bin_dir` (including `augustus_bin_dir`, `config_path` where applicable). Input files must exist; output paths must not.
- **Output:** A dictionary with `status`, `backend`, `version`, `command`, `return_code`, `outputs`, `log_tail`, `provenance` and parsed counts such as `lines`, `records`, `genes` or `transcripts`.
- **Example code:** not shown; the example needs data that is too large or an external tool that is not part of the AUGUSTUS package (see the tests and the manifest notes).

## 19) `dnakit.annotation.augustus_phastcons_to_hints()`

- **Purpose:** Convert a phastCons conserved-element table to hints (``phastconsDB2hints.pl``).
- **Method:** Runs `phastconsDB2hints.pl` through the DNAKit bounded executor.
- **API:** `dnakit.annotation.augustus_phastcons_to_hints(table_path, output_path, *, priority=None, source=None, min_score=None, cut_ends=None, augustus_bin_dir=None, verbose=True)`.
- **Input:** Required: `table_path`, `output_path`. Optional keyword-only: `priority`, `source`, `min_score`, `cut_ends`, `augustus_bin_dir` (including `augustus_bin_dir`, `config_path` where applicable). Input files must exist; output paths must not.
- **Output:** A dictionary with `status`, `backend`, `version`, `command`, `return_code`, `outputs`, `log_tail`, `provenance` and parsed counts such as `lines`, `records`, `genes` or `transcripts`.
- **Example code:** not shown; the example needs data that is too large or an external tool that is not part of the AUGUSTUS package (see the tests and the manifest notes).

## 20) `dnakit.annotation.augustus_polya_to_hints()`

- **Purpose:** Convert polyA/polyT positions to transcription-end (tts) hints.
- **Method:** Runs `hint_radius`, `new`, `newer`, `pasapolyA2hints.pl` through the DNAKit bounded executor. ``input_format="table"`` reads a polyA/polyT table (``polyA2hints.pl``; ``table_format`` ``new`` or ``newer``, ``hint_radius``, ``swap_strand`` for polyT files); ``input_format="pasa"`` reads PASA polyA output in FASTA format (``pasapolyA2hints.pl``, ``radius``).
- **API:** `dnakit.annotation.augustus_polya_to_hints(input_path, output_path, *, input_format='table', table_format=None, hint_radius=None, swap_strand=False, radius=None, augustus_bin_dir=None, verbose=True)`.
- **Input:** Required: `input_path`, `output_path`. Optional keyword-only: `input_format`, `table_format`, `hint_radius`, `swap_strand`, `radius`, `augustus_bin_dir` (including `augustus_bin_dir`, `config_path` where applicable). Input files must exist; output paths must not.
- **Output:** A dictionary with `status`, `backend`, `version`, `command`, `return_code`, `outputs`, `log_tail`, `provenance` and parsed counts such as `lines`, `records`, `genes` or `transcripts`.
- **Example code:** not shown; the example needs data that is too large or an external tool that is not part of the AUGUSTUS package (see the tests and the manifest notes).

## 21) `dnakit.annotation.augustus_psl_map()`

- **Purpose:** Chain PSL alignments transitively through a second PSL map (``pslMap.pl --in --map``).
- **Method:** Runs the AUGUSTUS program through the DNAKit bounded executor.
- **API:** `dnakit.annotation.augustus_psl_map(psl_path, map_path, output_path, *, augustus_bin_dir=None, verbose=True)`.
- **Input:** Required: `psl_path`, `map_path`, `output_path`. Optional keyword-only: `augustus_bin_dir` (including `augustus_bin_dir`, `config_path` where applicable). Input files must exist; output paths must not.
- **Output:** A dictionary with `status`, `backend`, `version`, `command`, `return_code`, `outputs`, `log_tail`, `provenance` and parsed counts such as `lines`, `records`, `genes` or `transcripts`.
- **Example code:** not shown; the example needs data that is too large or an external tool that is not part of the AUGUSTUS package (see the tests and the manifest notes).

## 22) `dnakit.annotation.augustus_psl_to_hints()`

- **Purpose:** Convert BLAT or GMAP PSL cDNA alignments to AUGUSTUS hints (``blat2hints.pl``).
- **Method:** Runs `blat2hints.pl` through the DNAKit bounded executor. The PSL must be sorted by target name (and start). ``column_offset=1`` reads the PSL layout of the UCSC database. Unset options keep the script defaults (priority 4, max gap 14, min intron 41, max intron 350000, source E).
- **API:** `dnakit.annotation.augustus_psl_to_hints(psl_path, output_path, *, priority=None, max_gap_length=None, min_intron_length=None, max_intron_length=None, max_query_gap_length=None, ep_cutoff=None, source=None, introns_only=False, no_multiplicity=False, remove_redundant=False, max_coverage=None, splice_site_hints=False, column_offset=None, score=None, augustus_bin_dir=None, verbose=True)`.
- **Input:** Required: `psl_path`, `output_path`. Optional keyword-only: `priority`, `max_gap_length`, `min_intron_length`, `max_intron_length`, `max_query_gap_length`, `ep_cutoff`, `source`, `introns_only`, `no_multiplicity`, `remove_redundant`, `max_coverage`, `splice_site_hints`, `column_offset`, `score`, `augustus_bin_dir` (including `augustus_bin_dir`, `config_path` where applicable). Input files must exist; output paths must not.
- **Output:** A dictionary with `status`, `backend`, `version`, `command`, `return_code`, `outputs`, `log_tail`, `provenance` and parsed counts such as `lines`, `records`, `genes` or `transcripts`.
- **Example code:**

```python
import random, tempfile
from pathlib import Path
from dnakit.annotation import augustus_p04_hints_functions as m
w = Path(tempfile.mkdtemp())
(w / "a.psl").write_text("100\t0\t0\t0\t0\t0\t1\t300\t+\tq1\t100\t0\t100\tchr1\t4000\t1000\t1400\t2\t50,50,\t0,50,\t1000,1350,\n")
r = m.augustus_psl_to_hints(str(w / "a.psl"), str(w / "psl.hints"), verbose=False)
print(r["lines"])
print(r["preview"][2])
```

- **Example output:**

```text
3
chr1	b2h	intron	1051	1350	0	.	.	grp=q1;pri=4;src=E
```

## 23) `dnakit.annotation.augustus_remove_redundant_hints()`

- **Purpose:** Remove hints from ``hints_path`` that are redundant with ``reference_hints_path`` (``rmRedundantHints.pl``).
- **Method:** Runs `hints_path`, `reference_hints_path`, `rmRedundantHints.pl` through the DNAKit bounded executor.
- **API:** `dnakit.annotation.augustus_remove_redundant_hints(hints_path, reference_hints_path, output_path, *, augustus_bin_dir=None, verbose=True)`.
- **Input:** Required: `hints_path`, `reference_hints_path`, `output_path`. Optional keyword-only: `augustus_bin_dir` (including `augustus_bin_dir`, `config_path` where applicable). Input files must exist; output paths must not.
- **Output:** A dictionary with `status`, `backend`, `version`, `command`, `return_code`, `outputs`, `log_tail`, `provenance` and parsed counts such as `lines`, `records`, `genes` or `transcripts`.
- **Example code:**

```python
import random, tempfile
from pathlib import Path
from dnakit.annotation import augustus_p04_hints_functions as m
w = Path(tempfile.mkdtemp())
(w / "h.gff").write_text("chr1\tb2h\tintron\t1051\t1350\t0\t+\t.\tgrp=a;pri=4;src=E\nchr1\tb2h\tintron\t1051\t1350\t0\t+\t.\tgrp=b;pri=4;src=E\n")
r = m.augustus_remove_redundant_hints(str(w / "h.gff"), str(w / "h.gff"), str(w / "nr.gff"), verbose=False)
print(r["status"], r["lines"])
```

- **Example output:**

```text
completed 1
```

## 24) `dnakit.annotation.augustus_retrogenes_to_hints()`

- **Purpose:** Convert a UCSC retroGenes table to hints (``retroDB2hints.pl``).
- **Method:** Runs `retroDB2hints.pl` through the DNAKit bounded executor.
- **API:** `dnakit.annotation.augustus_retrogenes_to_hints(table_path, output_path, *, priority=None, source=None, min_score=None, max_small_gap_length=None, augustus_bin_dir=None, verbose=True)`.
- **Input:** Required: `table_path`, `output_path`. Optional keyword-only: `priority`, `source`, `min_score`, `max_small_gap_length`, `augustus_bin_dir` (including `augustus_bin_dir`, `config_path` where applicable). Input files must exist; output paths must not.
- **Output:** A dictionary with `status`, `backend`, `version`, `command`, `return_code`, `outputs`, `log_tail`, `provenance` and parsed counts such as `lines`, `records`, `genes` or `transcripts`.
- **Example code:** not shown; the example needs data that is too large or an external tool that is not part of the AUGUSTUS package (see the tests and the manifest notes).

## 25) `dnakit.annotation.augustus_split_hal_to_maf()`

- **Purpose:** Export a HAL genome alignment to MAF chunks of a given size with overlap (``hal2maf_split.pl``).
- **Method:** Runs `hal2maf`, `hal2maf_split.pl`, `halStats`, `hal_bin_dir` through the DNAKit bounded executor. Needs the HAL tools (``hal2maf``, ``halStats``) in ``hal_bin_dir`` or on PATH. ``no_split_list_path`` lists genic intervals (seqname, start, end) inside which the alignment may not be split. The MAF files are written to the new ``output_dir``.
- **API:** `dnakit.annotation.augustus_split_hal_to_maf(hal_path, reference_genome, output_dir, *, reference_sequence=None, chunk_size=None, overlap=None, keep_dupes=False, keep_ancestors=False, no_split_list_path=None, cpus=None, hal_bin_dir=None, augustus_bin_dir=None, verbose=True)`.
- **Input:** Required: `hal_path`, `reference_genome`, `output_dir`. Optional keyword-only: `reference_sequence`, `chunk_size`, `overlap`, `keep_dupes`, `keep_ancestors`, `no_split_list_path`, `cpus`, `hal_bin_dir`, `augustus_bin_dir` (including `augustus_bin_dir`, `config_path` where applicable). Input files must exist; output paths must not.
- **Output:** A dictionary with `status`, `backend`, `version`, `command`, `return_code`, `outputs`, `log_tail`, `provenance` and parsed counts such as `lines`, `records`, `genes` or `transcripts`.
- **Example code:** not shown; the example needs data that is too large or an external tool that is not part of the AUGUSTUS package (see the tests and the manifest notes).

## 26) `dnakit.annotation.augustus_transmap_to_hints()`

- **Purpose:** Convert transMap alignments (genePred-like table) to hints (``transMap2hints.pl``); ``keep_ids_path`` restricts the genes.
- **Method:** Runs `keep_ids_path`, `transMap2hints.pl` through the DNAKit bounded executor.
- **API:** `dnakit.annotation.augustus_transmap_to_hints(genepred_path, output_path, *, keep_ids_path=None, priority=None, ep_cutoff=None, ep_margin=None, ip_cutoff=None, utr_end_cutoff=None, min_intron_length=None, min_intron_length_utr=None, start_stop_radius=None, tss_tts_radius=None, source=None, augustus_bin_dir=None, verbose=True)`.
- **Input:** Required: `genepred_path`, `output_path`. Optional keyword-only: `keep_ids_path`, `priority`, `ep_cutoff`, `ep_margin`, `ip_cutoff`, `utr_end_cutoff`, `min_intron_length`, `min_intron_length_utr`, `start_stop_radius`, `tss_tts_radius`, `source`, `augustus_bin_dir` (including `augustus_bin_dir`, `config_path` where applicable). Input files must exist; output paths must not.
- **Output:** A dictionary with `status`, `backend`, `version`, `command`, `return_code`, `outputs`, `log_tail`, `provenance` and parsed counts such as `lines`, `records`, `genes` or `transcripts`.
- **Example code:** not shown; the example needs data that is too large or an external tool that is not part of the AUGUSTUS package (see the tests and the manifest notes).

## 27) `dnakit.annotation.augustus_weed_maf()`

- **Purpose:** Remove columns that contain only gaps from the blocks of a MAF alignment (``weedMaf.pl``).
- **Method:** Runs `weedMaf.pl` through the DNAKit bounded executor.
- **API:** `dnakit.annotation.augustus_weed_maf(maf_path, output_path, *, augustus_bin_dir=None, verbose=True)`.
- **Input:** Required: `maf_path`, `output_path`. Optional keyword-only: `augustus_bin_dir` (including `augustus_bin_dir`, `config_path` where applicable). Input files must exist; output paths must not.
- **Output:** A dictionary with `status`, `backend`, `version`, `command`, `return_code`, `outputs`, `log_tail`, `provenance` and parsed counts such as `lines`, `records`, `genes` or `transcripts`.
- **Example code:** not shown; the example needs data that is too large or an external tool that is not part of the AUGUSTUS package (see the tests and the manifest notes).

## 28) `dnakit.annotation.augustus_wiggle_to_hints()`

- **Purpose:** Convert a wiggle coverage (or conservation) track to exonpart/CDSpart hints (``wig2hints.pl``).
- **Method:** Runs `wig2hints.pl` through the DNAKit bounded executor.
- **API:** `dnakit.annotation.augustus_wiggle_to_hints(wiggle_path, output_path, *, width=None, radius=None, margin=None, min_score=None, min_threshold=None, hint_type=None, prune=None, priority=None, strand=None, source=None, augustus_bin_dir=None, verbose=True)`.
- **Input:** Required: `wiggle_path`, `output_path`. Optional keyword-only: `width`, `radius`, `margin`, `min_score`, `min_threshold`, `hint_type`, `prune`, `priority`, `strand`, `source`, `augustus_bin_dir` (including `augustus_bin_dir`, `config_path` where applicable). Input files must exist; output paths must not.
- **Output:** A dictionary with `status`, `backend`, `version`, `command`, `return_code`, `outputs`, `log_tail`, `provenance` and parsed counts such as `lines`, `records`, `genes` or `transcripts`.
- **Example code:**

```python
import random, tempfile
from pathlib import Path
from dnakit.annotation import augustus_p04_hints_functions as m
w = Path(tempfile.mkdtemp())
(w / "cov.wig").write_text("variableStep chrom=chr1\n" + "".join(f"{p} 5\n" for p in range(101, 260)))
r = m.augustus_wiggle_to_hints(str(w / "cov.wig"), str(w / "wig.hints"), priority=2, verbose=False)
print(r["lines"] > 0, r["preview"][0].split("\t")[2])
```

- **Example output:**

```text
True CDSpart
```

