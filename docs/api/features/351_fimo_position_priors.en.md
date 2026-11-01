# FIMO motif scan with position-specific priors

This group has one function in `dnakit.regulation`: a typed command-line facade over `fimo` of the [MEME suite](https://meme-suite.org) (tested with 5.5.9). It scans DNA sequences for occurrences of MEME-format motifs and, when a position-specific prior (PSP) is given, scores each occurrence with the log-posterior odds of Cuellar-Partida et al. 2012 (Bioinformatics 28:56-62) instead of the plain log-odds. The prior can come from an epigenomic signal such as DNase or histone-mark accessibility, so a motif match in an open region outranks the same match in a closed region.

DNAKit does not bundle the MEME suite. Install it (`conda install -c bioconda meme`, or a release from meme-suite.org) and put `fimo` on `PATH`, or pass `fimo_path`. The prior distribution file is made by `compute-prior-dist` or `create-priors` of the same suite. If `fimo` is missing the function raises an error that includes the install hint. No shell is used and no network access happens.

The function covers the work-list operation "motif occurrence scan with position-specific priors" (`--psp` and `--prior-dist`). Plain PWM scanning without priors is already available as `scan_motifs_with_pwms()` in `dnakit.predictions`.

## Example data

A 60 bp sequence with the site `GATCAG` at positions 13, 22 and 34, a motif for `GATCAG`, a wiggle prior that raises the prior at position 13 only (positions without a value get the median prior), and the binned prior distribution (minimum, maximum, median, then the share of priors in each of four bins).

**Example file:** `example_seq.fa`

```text
>seq1
ACGTTGCAACCTGATCAGTAGGATCAGCTTAACGATCAGGTCCATTGACAGTCCGATGCA
```

**Example file:** `example_motif.meme`

```text
MEME version 4

ALPHABET= ACGT

strands: + -

Background letter frequencies
A 0.25 C 0.25 G 0.25 T 0.25

MOTIF TESTMOTIF
letter-probability matrix: alength= 4 w= 6 nsites= 20 E= 0
0.05 0.05 0.85 0.05
0.85 0.05 0.05 0.05
0.05 0.05 0.05 0.85
0.05 0.85 0.05 0.05
0.85 0.05 0.05 0.05
0.05 0.05 0.85 0.05
```

**Example file:** `example.wig`

```text
variableStep chrom=seq1 span=1
12 0.05
13 0.05
```

**Example file:** `example.dist`

```text
#min 0.0100
#max 0.0500
#median 0.0100
0.8
0.0
0.0
0.2
```

## 1) `dnakit.regulation.scan_motifs_with_position_priors()`

- **Purpose:** Find motif occurrences in DNA with FIMO, scored with position-specific priors (log-posterior odds), with p-values and q-values for that score.
- **Method:** Run `fimo --psp <prior_file> --prior-dist <prior_distribution_file> --alpha <prior_alpha> --thresh <p_value_threshold> ... <motif_file> <sequence_file>` into a new output folder, then parse `fimo.tsv`. Both prior files are required together, as FIMO requires. The defaults are those of FIMO: p-value threshold 1e-4, both strands, motif pseudocount 0.1, prior alpha 1.0, 100,000 stored scores. Options: `use_q_value_threshold` (`--qv-thresh`), `compute_q_values=False` (`--no-qvalue`), `scan_reverse_complement=False` (`--norc`), `max_strand`, `motif_ids` (`--motif`, repeated), `motif_pseudocount`, `background_file` or `background_source` (`motif`, `uniform` or `nrdb`, for `--bfile`), `parse_genomic_coordinates=False` (`--no-pgc`) and `max_stored_scores`. The wiggle coordinate of a prior is the 1-based start of the site.
- **API:** `dnakit.regulation.scan_motifs_with_position_priors(motif_file, sequence_file, prior_file, prior_distribution_file, output_dir, *, prior_alpha=1.0, p_value_threshold=1e-4, use_q_value_threshold=False, compute_q_values=True, scan_reverse_complement=True, max_strand=False, motif_ids=None, motif_pseudocount=0.1, background_file=None, background_source=None, parse_genomic_coordinates=True, max_stored_scores=100000, max_rows=1000, fimo_path=None)`.
- **Input:** A MEME-format motif file, a FASTA file, a prior in wiggle or MEME PSP format, its binned distribution, and a folder that does not exist yet.
- **Output:** A dictionary with `status` (`completed` or `failed`), `backend`, `version`, `command`, `return_code`, `log_tail`, `output_dir`, `outputs` (`fimo.tsv`, `fimo.gff`, `fimo.html`, `fimo.xml`, `cisml.xml`, `best_site.narrowPeak`), `match_count`, `matches` (the first `max_rows` rows with `motif_id`, `motif_alt_id`, `sequence_name`, `start`, `stop`, `strand`, `score`, `p_value`, `q_value`, `matched_sequence`) and `matches_truncated`.
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.regulation.fimo_p17_cli_functions import scan_motifs_with_position_priors

output = Path(tempfile.mkdtemp()) / "fimo_out"
result = scan_motifs_with_position_priors(
    "example_motif.meme",
    "example_seq.fa",
    "example.wig",
    "example.dist",
    str(output),
    p_value_threshold=1e-3,
    verbose=False,
)
print(result["status"], result["backend"], result["match_count"])
for row in result["matches"]:
    print(row["motif_id"], row["start"], row["stop"], row["strand"], row["score"], row["p_value"], row["matched_sequence"])
print(sorted(Path(path).name for path in result["outputs"]))
```

- **Example output:**

```text
completed fimo 4
TESTMOTIF 13 18 + 6.35714 4.74e-05 GATCAG
TESTMOTIF 11 16 - 3.98214 0.000237 GATCAG
TESTMOTIF 22 27 + 3.98214 0.000237 GATCAG
TESTMOTIF 34 39 + 3.98214 0.000237 GATCAG
['best_site.narrowPeak', 'cisml.xml', 'fimo.gff', 'fimo.html', 'fimo.tsv', 'fimo.xml']
```
