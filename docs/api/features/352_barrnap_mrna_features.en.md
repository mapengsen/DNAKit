# Barrnap messenger-RNA feature prediction

This group has one function in `dnakit.annotation`: a typed command-line facade over `barrnap --mrna` of [Barrnap](https://github.com/tseemann/barrnap) 1.x (written against 1.10.5). Besides rRNA, Barrnap 1.x predicts protein-coding genes with Pyrodigal and reports, as GFF3, the `mRNA` (extended upstream over the ribosome binding site), the `RBS` (with its motif), the `CDS` (named by a DIAMOND search against the Barrnap protein database of the chosen kingdom), the `sig_pep` (signal peptide, ZiggyPep) and the Rho-independent `terminator` (TransTermHP).

DNAKit does not bundle Barrnap. Install it with its dependencies (`conda install -c bioconda -c conda-forge barrnap`) and put `barrnap` on `PATH`, or pass `barrnap_path`. If it is missing the function raises an error that includes the install hint. No shell is used and no network access happens.

Barrnap skips the mRNA search with a warning when the combined input is shorter than about 19 kb, so give it genomes or contigs, not single genes. The rRNA, tRNA and ncRNA searches of Barrnap are separate operations; `include_rrna=True` only adds Barrnap's default rRNA search to the same GFF3 file.

## Example data

Any bacterial genome or contig FASTA of at least 19 kb, for example `genome.fa`. The example below does not declare it, so the page is not run by the documentation checker; the output shows the shape of a real result, not measured values.

## 1) `dnakit.annotation.predict_mrna_features_barrnap()`

- **Purpose:** Predict protein-coding genes with RBS, signal peptide and terminator features in DNA with Barrnap.
- **Method:** Run `barrnap --quiet --mrna --no-rrna --kingdom <kingdom> --gcode <genetic_code> --evalue <evalue> --threads <threads> <sequence files>`, stream the GFF3 from standard output into the new `output_gff`, keep standard error separate and parse the features. `kingdom` is `bac`, `arc` or `fun`; with `genetic_code='auto'` Barrnap uses table 11, 25 or 1 respectively. `add_ids` is `--addids`, `include_rrna` replaces `--no-rrna` by `--rrna`, `database_dir` is `--dbdir`.
- **API:** `dnakit.annotation.predict_mrna_features_barrnap(sequence_files, output_gff, *, kingdom="bac", genetic_code="auto", evalue=0.001, threads=1, add_ids=False, include_rrna=False, database_dir=None, max_features=1000, barrnap_path=None)`.
- **Input:** One or more DNA files (FASTA, GenBank or FASTQ, optionally gzip-compressed) and a GFF3 path that does not exist yet.
- **Output:** A dictionary with `status`, `backend`, `version`, `command`, `return_code`, `log_tail`, `gff_file`, `feature_counts` (type to count), `feature_count`, `features` (the first `max_features` rows with `seqid`, `source`, `type`, `start`, `end`, `score`, `strand`, `phase` and an `attributes` dictionary) and `features_truncated`.
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.annotation.barrnap_p17_cli_functions import predict_mrna_features_barrnap

output = Path(tempfile.mkdtemp()) / "mrna.gff3"
result = predict_mrna_features_barrnap(["genome.fa"], str(output), kingdom="bac", threads=4)
print(result["status"], result["feature_counts"])
```

- **Example output (shape only, not measured):**

```text
completed {'mRNA': 190, 'RBS': 170, 'CDS': 190, 'sig_pep': 12, 'terminator': 40}
```
