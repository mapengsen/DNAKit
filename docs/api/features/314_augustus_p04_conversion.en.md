# AUGUSTUS: gene-model, sequence and masking conversions

The functions in `dnakit.annotation` listed here are typed facades over AUGUSTUS 3.5.0 programs and helper scripts ([AUGUSTUS](https://github.com/Gaius-Augustus/Augustus)). Each one runs a single fixed program with an explicit argument list (no shell). Install AUGUSTUS (`conda install -c bioconda augustus`) and put its `bin` folder on PATH, or pass `augustus_bin_dir`. If the program is missing the function raises `BackendUnavailableError` with the install hint. The species parameter sets are read from `AUGUSTUS_CONFIG_PATH` (or `config_path`); training functions write there, so the folder must be writable. Output paths must be new; nothing is overwritten. Each result is a dictionary with `status`, `backend`, `version`, `command`, `return_code`, `outputs`, `log_tail` and `provenance`, plus small parsed counts.

## 1) `dnakit.annotation.augustus_add_names_to_gff3()`

- **Purpose:** Add ``Name=`` attributes to the gene and mRNA lines of an AUGUSTUS GFF3 (``add_name_to_gff3.pl``); ``only_gene_structure`` keeps gene, mRNA, CDS and exon only.
- **Method:** Runs `add_name_to_gff3.pl`, `only_gene_structure` through the DNAKit bounded executor.
- **API:** `dnakit.annotation.augustus_add_names_to_gff3(gff3_path, output_path, *, only_gene_structure=False, augustus_bin_dir=None, verbose=True)`.
- **Input:** Required: `gff3_path`, `output_path`. Optional keyword-only: `only_gene_structure`, `augustus_bin_dir` (including `augustus_bin_dir`, `config_path` where applicable). Input files must exist; output paths must not.
- **Output:** A dictionary with `status`, `backend`, `version`, `command`, `return_code`, `outputs`, `log_tail`, `provenance` and parsed counts such as `lines`, `records`, `genes` or `transcripts`.
- **Example code:**

```python
import random, tempfile
from pathlib import Path
from dnakit.annotation import augustus_p04_conversion_functions as m
w = Path(tempfile.mkdtemp())
(w / "genes.gtf").write_text(
    'chr1\tAUGUSTUS\tgene\t101\t400\t.\t+\t.\tg1\n'
    'chr1\tAUGUSTUS\ttranscript\t101\t400\t.\t+\t.\tg1.t1\n'
    'chr1\tAUGUSTUS\tstart_codon\t101\t103\t.\t+\t0\ttranscript_id "g1.t1"; gene_id "g1";\n'
    'chr1\tAUGUSTUS\tCDS\t101\t200\t.\t+\t0\ttranscript_id "g1.t1"; gene_id "g1";\n'
    'chr1\tAUGUSTUS\tCDS\t301\t400\t.\t+\t1\ttranscript_id "g1.t1"; gene_id "g1";\n'
)
g = m.augustus_gtf_to_gff(str(w / "genes.gtf"), str(w / "genes.gff3"), gff3=True, print_exons=True, verbose=False)
r = m.augustus_add_names_to_gff3(g["outputs"]["output"], str(w / "named.gff3"), verbose=False)
print(r["preview"][1])
```

- **Example output:**

```text
chr1	AUGUSTUS	mRNA	101	400	.	+	.	ID=g1.t1;Parent=g1;Name=g1.t1
```

## 2) `dnakit.annotation.augustus_cegma_to_gff()`

- **Purpose:** Convert a CEGMA output GFF to the CDS GTF used for AUGUSTUS training (``cegma2gff.pl``).
- **Method:** Runs `cegma2gff.pl` through the DNAKit bounded executor.
- **API:** `dnakit.annotation.augustus_cegma_to_gff(cegma_gff_path, output_path, *, augustus_bin_dir=None, verbose=True)`.
- **Input:** Required: `cegma_gff_path`, `output_path`. Optional keyword-only: `augustus_bin_dir` (including `augustus_bin_dir`, `config_path` where applicable). Input files must exist; output paths must not.
- **Output:** A dictionary with `status`, `backend`, `version`, `command`, `return_code`, `outputs`, `log_tail`, `provenance` and parsed counts such as `lines`, `records`, `genes` or `transcripts`.
- **Example code:**

```python
import random, tempfile
from pathlib import Path
from dnakit.annotation import augustus_p04_conversion_functions as m
w = Path(tempfile.mkdtemp())
(w / "cegma.gff").write_text("CHR\tcegma\tFirst\t100\t200\t5.0\t+\t0\tK1.1\nCHR\tcegma\tExon\t100\t200\t5.0\t+\t.\tK1.1\n")
r = m.augustus_cegma_to_gff(str(w / "cegma.gff"), str(w / "cds.gtf"), verbose=False)
print(r["preview"][0])
```

- **Example output:**

```text
CHR	cegma	CDS	100	200	5.0	+	0	transcript_id "g1.K1.1"
```

## 3) `dnakit.annotation.augustus_check_gff_names_in_fasta()`

- **Purpose:** Check that the sequence names of a GFF are FASTA headers of the genome and that it has 9 columns (``findGffNamesInFasta.pl``).
- **Method:** Runs `findGffNamesInFasta.pl`, `names_missing`, `output_path`, `wrong_column_count` through the DNAKit bounded executor. ``names_missing`` is True when the GFF uses sequence names absent from the genome, ``wrong_column_count`` is True when a line does not have 9 columns (the script writes this flag to ``output_path``).
- **API:** `dnakit.annotation.augustus_check_gff_names_in_fasta(gff_path, genome_path, output_path, *, augustus_bin_dir=None, verbose=True)`.
- **Input:** Required: `gff_path`, `genome_path`, `output_path`. Optional keyword-only: `augustus_bin_dir` (including `augustus_bin_dir`, `config_path` where applicable). Input files must exist; output paths must not.
- **Output:** A dictionary with `status`, `backend`, `version`, `command`, `return_code`, `outputs`, `log_tail`, `provenance` and parsed counts such as `lines`, `records`, `genes` or `transcripts`.
- **Example code:**

```python
import random, tempfile
from pathlib import Path
from dnakit.annotation import augustus_p04_conversion_functions as m
w = Path(tempfile.mkdtemp())
rng = random.Random(7)
(w / "genome.fa").write_text(">chr1\n" + "".join(rng.choices("ACGT", k=4000)) + "\n")
(w / "x.gff").write_text("chr1\t.\tgene\t1\t2\t.\t+\t.\tID=a\nchrX\t.\tgene\t1\t2\t.\t+\t.\tID=b\n")
r = m.augustus_check_gff_names_in_fasta(str(w / "x.gff"), str(w / "genome.fa"), str(w / "chk.txt"), verbose=False)
print(r["names_missing"], r["wrong_column_count"])
```

- **Example output:**

```text
True False
```

## 4) `dnakit.annotation.augustus_cluster_gene_dependencies()`

- **Purpose:** Group genes into clusters without inter-cluster dependencies from a pair-dependency file (``clusterAndSplitGenes.pl``); clusters go to ``clusters_path``.
- **Method:** Runs `clusterAndSplitGenes.pl`, `clusters_path` through the DNAKit bounded executor.
- **API:** `dnakit.annotation.augustus_cluster_gene_dependencies(dependency_path, clusters_path, *, threshold=None, max_cluster_size=None, transcript_format=False, augustus_bin_dir=None, verbose=True)`.
- **Input:** Required: `dependency_path`, `clusters_path`. Optional keyword-only: `threshold`, `max_cluster_size`, `transcript_format`, `augustus_bin_dir` (including `augustus_bin_dir`, `config_path` where applicable). Input files must exist; output paths must not.
- **Output:** A dictionary with `status`, `backend`, `version`, `command`, `return_code`, `outputs`, `log_tail`, `provenance` and parsed counts such as `lines`, `records`, `genes` or `transcripts`.
- **Example code:** not shown; the example needs data that is too large or an external tool that is not part of the AUGUSTUS package (see the tests and the manifest notes).

## 5) `dnakit.annotation.augustus_compare_masking()`

- **Purpose:** Compare the soft-masking of two versions of the same assembly (``compare_masking.pl``; both are held in memory). The report goes to ``output_path``.
- **Method:** Runs `compare_masking.pl`, `output_path` through the DNAKit bounded executor.
- **API:** `dnakit.annotation.augustus_compare_masking(softmasked_path, second_softmasked_path, output_path, *, augustus_bin_dir=None, verbose=True)`.
- **Input:** Required: `softmasked_path`, `second_softmasked_path`, `output_path`. Optional keyword-only: `augustus_bin_dir` (including `augustus_bin_dir`, `config_path` where applicable). Input files must exist; output paths must not.
- **Output:** A dictionary with `status`, `backend`, `version`, `command`, `return_code`, `outputs`, `log_tail`, `provenance` and parsed counts such as `lines`, `records`, `genes` or `transcripts`.
- **Example code:**

```python
import random, tempfile
from pathlib import Path
from dnakit.annotation import augustus_p04_conversion_functions as m
w = Path(tempfile.mkdtemp())
rng = random.Random(7)
(w / "genome.fa").write_text(">chr1\n" + "".join(rng.choices("ACGT", k=4000)) + "\n")
r = m.augustus_compare_masking(str(w / "genome.fa"), str(w / "genome.fa"), str(w / "cmp.txt"), verbose=False)
print(r["preview"])
```

- **Example output:**

```text
['Masked in both files: 0', 'Masked in File1, only: 0', 'Masked in File2, only: 0']
```

## 6) `dnakit.annotation.augustus_compute_flanking_region()`

- **Purpose:** Average gene length and the flanking region derived from a GTF (``computeFlankingRegion.pl``); the script text is written to ``output_path`` and the numbers are parsed into ``numbers``.
- **Method:** Runs `computeFlankingRegion.pl`, `numbers`, `output_path` through the DNAKit bounded executor.
- **API:** `dnakit.annotation.augustus_compute_flanking_region(gtf_path, output_path, *, augustus_bin_dir=None, verbose=True)`.
- **Input:** Required: `gtf_path`, `output_path`. Optional keyword-only: `augustus_bin_dir` (including `augustus_bin_dir`, `config_path` where applicable). Input files must exist; output paths must not.
- **Output:** A dictionary with `status`, `backend`, `version`, `command`, `return_code`, `outputs`, `log_tail`, `provenance` and parsed counts such as `lines`, `records`, `genes` or `transcripts`.
- **Example code:**

```python
import random, tempfile
from pathlib import Path
from dnakit.annotation import augustus_p04_conversion_functions as m
w = Path(tempfile.mkdtemp())
(w / "genes.gtf").write_text(
    'chr1\tAUGUSTUS\tgene\t101\t400\t.\t+\t.\tg1\n'
    'chr1\tAUGUSTUS\ttranscript\t101\t400\t.\t+\t.\tg1.t1\n'
    'chr1\tAUGUSTUS\tstart_codon\t101\t103\t.\t+\t0\ttranscript_id "g1.t1"; gene_id "g1";\n'
    'chr1\tAUGUSTUS\tCDS\t101\t200\t.\t+\t0\ttranscript_id "g1.t1"; gene_id "g1";\n'
    'chr1\tAUGUSTUS\tCDS\t301\t400\t.\t+\t1\ttranscript_id "g1.t1"; gene_id "g1";\n'
)
r = m.augustus_compute_flanking_region(str(w / "genes.gtf"), str(w / "flank.txt"), verbose=False)
print(r["preview"][1])
```

- **Example output:**

```text
The flanking_DNA value is: 150 (the Minimum of 10 000 and 150)
```

## 7) `dnakit.annotation.augustus_convert_genepred()`

- **Purpose:** Convert UCSC ``genePredExt`` rows to the bigGenePred-ready layout (``gp2othergp.pl``).
- **Method:** Runs `genePredExt`, `gp2othergp.pl` through the DNAKit bounded executor.
- **API:** `dnakit.annotation.augustus_convert_genepred(genepred_path, output_path, *, augustus_bin_dir=None, verbose=True)`.
- **Input:** Required: `genepred_path`, `output_path`. Optional keyword-only: `augustus_bin_dir` (including `augustus_bin_dir`, `config_path` where applicable). Input files must exist; output paths must not.
- **Output:** A dictionary with `status`, `backend`, `version`, `command`, `return_code`, `outputs`, `log_tail`, `provenance` and parsed counts such as `lines`, `records`, `genes` or `transcripts`.
- **Example code:**

```python
import random, tempfile
from pathlib import Path
from dnakit.annotation import augustus_p04_conversion_functions as m
w = Path(tempfile.mkdtemp())
(w / "gp.txt").write_text("g1\tchr1\t+\t100\t900\t0\t1\t2\t100,500\t300,900\n")
r = m.augustus_convert_genepred(str(w / "gp.txt"), str(w / "big.gp"), verbose=False)
print(r["preview"][0])
```

- **Example output:**

```text
chr1	100	900	g1	0	+	0	1	255,0,0	2	200,400	0,400	g1				coding	?	?	?
```

## 8) `dnakit.annotation.augustus_extract_annotation_fasta()`

- **Purpose:** Write protein and coding sequences from AUGUSTUS output comments (``getAnnoFasta.pl``).
- **Method:** Runs `.cdsexons`, `.codingseq`, `.mrna`, `augustus.aa` through the DNAKit bounded executor. The prediction must contain ``# protein sequence`` and ``# coding sequence`` blocks (``protein_sequences=True`` and ``coding_sequences=True`` in ``augustus_predict_genes``). With ``genome_path`` the exon sequences (``.cdsexons``) and the mRNA including UTRs (``.mrna``) are written too. Files are ``augustus.aa``, ``.codingseq`` and so on in the new ``output_dir``.
- **API:** `dnakit.annotation.augustus_extract_annotation_fasta(prediction_path, output_dir, *, genome_path=None, chop_cds=False, augustus_bin_dir=None, verbose=True)`.
- **Input:** Required: `prediction_path`, `output_dir`. Optional keyword-only: `genome_path`, `chop_cds`, `augustus_bin_dir` (including `augustus_bin_dir`, `config_path` where applicable). Input files must exist; output paths must not.
- **Output:** A dictionary with `status`, `backend`, `version`, `command`, `return_code`, `outputs`, `log_tail`, `provenance` and parsed counts such as `lines`, `records`, `genes` or `transcripts`.
- **Example code:** not shown; the example needs data that is too large or an external tool that is not part of the AUGUSTUS package (see the tests and the manifest notes).

## 9) `dnakit.annotation.augustus_extract_cds_sequences()`

- **Purpose:** Extract coding sequences and translated proteins of CDS models in a GTF/GFF3 plus genome FASTA.
- **Method:** Runs `drop_invalid_stops`, `getAnnoFastaFromJoingenes.py`, `gtf2aa.pl`, `stringtie2fa.py` through the DNAKit bounded executor. ``method="getAnnoFastaFromJoingenes"`` (``getAnnoFastaFromJoingenes.py``) writes ``<stem>.codingseq`` and ``<stem>.aa`` (GTF or GFF3, ``translation_table``, ``drop_invalid_stops``); ``"stringtie2fa"`` (``stringtie2fa.py``, GTF only) is the earlier variant with the same outputs; ``"gtf2aa"`` (``gtf2aa.pl``) writes only the protein FASTA ``<stem>.aa``.
- **API:** `dnakit.annotation.augustus_extract_cds_sequences(genome_path, annotation_path, output_stem, *, method='getAnnoFastaFromJoingenes', annotation_format='gtf', translation_table=None, drop_invalid_stops=False, augustus_bin_dir=None, verbose=True)`.
- **Input:** Required: `genome_path`, `annotation_path`, `output_stem`. Optional keyword-only: `method`, `annotation_format`, `translation_table`, `drop_invalid_stops`, `augustus_bin_dir` (including `augustus_bin_dir`, `config_path` where applicable). Input files must exist; output paths must not.
- **Output:** A dictionary with `status`, `backend`, `version`, `command`, `return_code`, `outputs`, `log_tail`, `provenance` and parsed counts such as `lines`, `records`, `genes` or `transcripts`.
- **Example code:**

```python
import random, tempfile
from pathlib import Path
from dnakit.annotation import augustus_p04_conversion_functions as m
w = Path(tempfile.mkdtemp())
rng = random.Random(7)
(w / "genome.fa").write_text(">chr1\n" + "".join(rng.choices("ACGT", k=4000)) + "\n")
(w / "genes.gtf").write_text(
    'chr1\tAUGUSTUS\tgene\t101\t400\t.\t+\t.\tg1\n'
    'chr1\tAUGUSTUS\ttranscript\t101\t400\t.\t+\t.\tg1.t1\n'
    'chr1\tAUGUSTUS\tstart_codon\t101\t103\t.\t+\t0\ttranscript_id "g1.t1"; gene_id "g1";\n'
    'chr1\tAUGUSTUS\tCDS\t101\t200\t.\t+\t0\ttranscript_id "g1.t1"; gene_id "g1";\n'
    'chr1\tAUGUSTUS\tCDS\t301\t400\t.\t+\t1\ttranscript_id "g1.t1"; gene_id "g1";\n'
)
r = m.augustus_extract_cds_sequences(str(w / "genome.fa"), str(w / "genes.gtf"), str(w / "prot"), method="gtf2aa", verbose=False)
print(r["records"])
```

- **Example output:**

```text
{'protein': 1}
```

## 10) `dnakit.annotation.augustus_extract_mrna_sequences()`

- **Purpose:** Splice mRNA sequences of GFF gene models out of the genome (``gffGetmRNA.pl``) into the FASTA ``output_path``.
- **Method:** Runs `gffGetmRNA.pl`, `output_path` through the DNAKit bounded executor.
- **API:** `dnakit.annotation.augustus_extract_mrna_sequences(gff_path, genome_path, output_path, *, augustus_bin_dir=None, verbose=True)`.
- **Input:** Required: `gff_path`, `genome_path`, `output_path`. Optional keyword-only: `augustus_bin_dir` (including `augustus_bin_dir`, `config_path` where applicable). Input files must exist; output paths must not.
- **Output:** A dictionary with `status`, `backend`, `version`, `command`, `return_code`, `outputs`, `log_tail`, `provenance` and parsed counts such as `lines`, `records`, `genes` or `transcripts`.
- **Example code:** not shown; the example needs data that is too large or an external tool that is not part of the AUGUSTUS package (see the tests and the manifest notes).

## 11) `dnakit.annotation.augustus_filter_genbank_genes()`

- **Purpose:** Select gene models of a GenBank file by a gene-name list.
- **Method:** Runs `exclude_loci`, `exclude_mrna`, `filterGenes.pl`, `filterGenesIn.pl` through the DNAKit bounded executor. ``mode`` ``exclude_loci`` drops the loci named in ``names_path`` (``filterGenes.pl``), ``include_loci`` keeps only them, first of identical loci (``filterGenesIn.pl``), ``include_mrna_gtf`` keeps the genes of a GTF (``names_path`` is then a GTF, ``filterGenesIn_mRNAname.pl``) and ``exclude_mrna`` drops listed mRNA names (``filterGenesOut_mRNAname.pl``).
- **API:** `dnakit.annotation.augustus_filter_genbank_genes(genbank_path, names_path, output_path, *, mode='exclude_loci', augustus_bin_dir=None, verbose=True)`.
- **Input:** Required: `genbank_path`, `names_path`, `output_path`. Optional keyword-only: `mode`, `augustus_bin_dir` (including `augustus_bin_dir`, `config_path` where applicable). Input files must exist; output paths must not.
- **Output:** A dictionary with `status`, `backend`, `version`, `command`, `return_code`, `outputs`, `log_tail`, `provenance` and parsed counts such as `lines`, `records`, `genes` or `transcripts`.
- **Example code:** not shown; the example needs data that is too large or an external tool that is not part of the AUGUSTUS package (see the tests and the manifest notes).

## 12) `dnakit.annotation.augustus_fix_joingenes_gtf()`

- **Purpose:** Add the missing ``gene`` feature lines to a ``joingenes`` GTF (``fix_joingenes_gtf.pl``).
- **Method:** Runs `fix_joingenes_gtf.pl`, `gene`, `joingenes` through the DNAKit bounded executor.
- **API:** `dnakit.annotation.augustus_fix_joingenes_gtf(gtf_path, output_path, *, augustus_bin_dir=None, verbose=True)`.
- **Input:** Required: `gtf_path`, `output_path`. Optional keyword-only: `augustus_bin_dir` (including `augustus_bin_dir`, `config_path` where applicable). Input files must exist; output paths must not.
- **Output:** A dictionary with `status`, `backend`, `version`, `command`, `return_code`, `outputs`, `log_tail`, `provenance` and parsed counts such as `lines`, `records`, `genes` or `transcripts`.
- **Example code:**

```python
import random, tempfile
from pathlib import Path
from dnakit.annotation import augustus_p04_conversion_functions as m
w = Path(tempfile.mkdtemp())
(w / "genes.gtf").write_text(
    'chr1\tAUGUSTUS\tgene\t101\t400\t.\t+\t.\tg1\n'
    'chr1\tAUGUSTUS\ttranscript\t101\t400\t.\t+\t.\tg1.t1\n'
    'chr1\tAUGUSTUS\tstart_codon\t101\t103\t.\t+\t0\ttranscript_id "g1.t1"; gene_id "g1";\n'
    'chr1\tAUGUSTUS\tCDS\t101\t200\t.\t+\t0\ttranscript_id "g1.t1"; gene_id "g1";\n'
    'chr1\tAUGUSTUS\tCDS\t301\t400\t.\t+\t1\ttranscript_id "g1.t1"; gene_id "g1";\n'
)
j = __import__("dnakit.annotation.augustus_p04_prediction_functions", fromlist=["x"]).augustus_join_gene_sets([str(w / "genes.gtf")], str(w / "joined.gtf"), verbose=False)
r = m.augustus_fix_joingenes_gtf(j["outputs"]["merged"], str(w / "fixed.gtf"), verbose=False)
print(r["preview"][0].split("\t")[2])
```

- **Example output:**

```text
gene
```

## 13) `dnakit.annotation.augustus_gbrowse_old_to_gff3()`

- **Purpose:** Convert legacy GBrowse-format AUGUSTUS predictions to GFF3 (``gbrowseold2gff3.pl``); ``hints`` reads hints instead of genes.
- **Method:** Runs `gbrowseold2gff3.pl`, `hints` through the DNAKit bounded executor.
- **API:** `dnakit.annotation.augustus_gbrowse_old_to_gff3(prediction_path, output_path, *, hints=False, source=None, augustus_bin_dir=None, verbose=True)`.
- **Input:** Required: `prediction_path`, `output_path`. Optional keyword-only: `hints`, `source`, `augustus_bin_dir` (including `augustus_bin_dir`, `config_path` where applicable). Input files must exist; output paths must not.
- **Output:** A dictionary with `status`, `backend`, `version`, `command`, `return_code`, `outputs`, `log_tail`, `provenance` and parsed counts such as `lines`, `records`, `genes` or `transcripts`.
- **Example code:** not shown; the example needs data that is too large or an external tool that is not part of the AUGUSTUS package (see the tests and the manifest notes).

## 14) `dnakit.annotation.augustus_genbank_loci()`

- **Purpose:** Extract the tab-separated transcript-to-locus list of a GenBank training file (``get_loci_from_gb.pl``).
- **Method:** Runs `get_loci_from_gb.pl`, `output_path` through the DNAKit bounded executor. The script prints the list to standard output and leaves its ``--out`` file empty (upstream bug); this facade points ``--out`` at a scratch file and writes the printed list to the new ``output_path``.
- **API:** `dnakit.annotation.augustus_genbank_loci(genbank_path, output_path, *, augustus_bin_dir=None, verbose=True)`.
- **Input:** Required: `genbank_path`, `output_path`. Optional keyword-only: `augustus_bin_dir` (including `augustus_bin_dir`, `config_path` where applicable). Input files must exist; output paths must not.
- **Output:** A dictionary with `status`, `backend`, `version`, `command`, `return_code`, `outputs`, `log_tail`, `provenance` and parsed counts such as `lines`, `records`, `genes` or `transcripts`.
- **Example code:**

```python
import random, tempfile
from pathlib import Path
from dnakit.annotation import augustus_p04_conversion_functions as m
w = Path(tempfile.mkdtemp())
rng = random.Random(7)
(w / "genome.fa").write_text(">chr1\n" + "".join(rng.choices("ACGT", k=4000)) + "\n")
(w / "genes.gtf").write_text(
    'chr1\tAUGUSTUS\tgene\t101\t400\t.\t+\t.\tg1\n'
    'chr1\tAUGUSTUS\ttranscript\t101\t400\t.\t+\t.\tg1.t1\n'
    'chr1\tAUGUSTUS\tstart_codon\t101\t103\t.\t+\t0\ttranscript_id "g1.t1"; gene_id "g1";\n'
    'chr1\tAUGUSTUS\tCDS\t101\t200\t.\t+\t0\ttranscript_id "g1.t1"; gene_id "g1";\n'
    'chr1\tAUGUSTUS\tCDS\t301\t400\t.\t+\t1\ttranscript_id "g1.t1"; gene_id "g1";\n'
)
gb = m.augustus_gff_to_genbank(str(w / "genes.gtf"), str(w / "genome.fa"), 50, str(w / "train.gb"), verbose=False)
r = m.augustus_genbank_loci(gb["outputs"]["output"], str(w / "loci.tsv"), verbose=False)
print(r["preview"])
```

- **Example output:**

```text
['g1.t1\tchr1_51-450']
```

## 15) `dnakit.annotation.augustus_genbank_to_gff()`

- **Purpose:** Convert GenBank training genes (with flanking-region LOCUS names such as ``chr1_100-900``) to genome-scale GFF (``gbSmallDNA2gff.pl``).
- **Method:** Runs `gbSmallDNA2gff.pl` through the DNAKit bounded executor.
- **API:** `dnakit.annotation.augustus_genbank_to_gff(genbank_path, output_path, *, augustus_bin_dir=None, verbose=True)`.
- **Input:** Required: `genbank_path`, `output_path`. Optional keyword-only: `augustus_bin_dir` (including `augustus_bin_dir`, `config_path` where applicable). Input files must exist; output paths must not.
- **Output:** A dictionary with `status`, `backend`, `version`, `command`, `return_code`, `outputs`, `log_tail`, `provenance` and parsed counts such as `lines`, `records`, `genes` or `transcripts`.
- **Example code:**

```python
import random, tempfile
from pathlib import Path
from dnakit.annotation import augustus_p04_conversion_functions as m
w = Path(tempfile.mkdtemp())
rng = random.Random(7)
(w / "genome.fa").write_text(">chr1\n" + "".join(rng.choices("ACGT", k=4000)) + "\n")
(w / "genes.gtf").write_text(
    'chr1\tAUGUSTUS\tgene\t101\t400\t.\t+\t.\tg1\n'
    'chr1\tAUGUSTUS\ttranscript\t101\t400\t.\t+\t.\tg1.t1\n'
    'chr1\tAUGUSTUS\tstart_codon\t101\t103\t.\t+\t0\ttranscript_id "g1.t1"; gene_id "g1";\n'
    'chr1\tAUGUSTUS\tCDS\t101\t200\t.\t+\t0\ttranscript_id "g1.t1"; gene_id "g1";\n'
    'chr1\tAUGUSTUS\tCDS\t301\t400\t.\t+\t1\ttranscript_id "g1.t1"; gene_id "g1";\n'
)
gb = m.augustus_gff_to_genbank(str(w / "genes.gtf"), str(w / "genome.fa"), 50, str(w / "train.gb"), verbose=False)
r = m.augustus_genbank_to_gff(gb["outputs"]["output"], str(w / "back.gff"), verbose=False)
print(r["lines"], r["preview"][0])
```

- **Example output:**

```text
2 chr1	someSource	CDS	101	200	.	+	.	gene_id "g1.t1"; transcript_id "g1.t1";
```

## 16) `dnakit.annotation.augustus_gff_to_genbank()`

- **Purpose:** Convert GFF gene models and a genome FASTA to a minimal GenBank training file with flanking DNA (``gff2gbSmallDNA.pl``).
- **Method:** Runs `gff2gbSmallDNA.pl` through the DNAKit bounded executor.
- **API:** `dnakit.annotation.augustus_gff_to_genbank(gff_path, genome_path, flanking_dna, output_path, *, exclude_genes_path=None, include_genes_path=None, allow_overlaps=False, connected=False, keep_softmasking=False, hardmask=False, augustus_bin_dir=None, verbose=True)`.
- **Input:** Required: `gff_path`, `genome_path`, `flanking_dna`, `output_path`. Optional keyword-only: `exclude_genes_path`, `include_genes_path`, `allow_overlaps`, `connected`, `keep_softmasking`, `hardmask`, `augustus_bin_dir` (including `augustus_bin_dir`, `config_path` where applicable). Input files must exist; output paths must not.
- **Output:** A dictionary with `status`, `backend`, `version`, `command`, `return_code`, `outputs`, `log_tail`, `provenance` and parsed counts such as `lines`, `records`, `genes` or `transcripts`.
- **Example code:**

```python
import random, tempfile
from pathlib import Path
from dnakit.annotation import augustus_p04_conversion_functions as m
w = Path(tempfile.mkdtemp())
rng = random.Random(7)
(w / "genome.fa").write_text(">chr1\n" + "".join(rng.choices("ACGT", k=4000)) + "\n")
(w / "genes.gtf").write_text(
    'chr1\tAUGUSTUS\tgene\t101\t400\t.\t+\t.\tg1\n'
    'chr1\tAUGUSTUS\ttranscript\t101\t400\t.\t+\t.\tg1.t1\n'
    'chr1\tAUGUSTUS\tstart_codon\t101\t103\t.\t+\t0\ttranscript_id "g1.t1"; gene_id "g1";\n'
    'chr1\tAUGUSTUS\tCDS\t101\t200\t.\t+\t0\ttranscript_id "g1.t1"; gene_id "g1";\n'
    'chr1\tAUGUSTUS\tCDS\t301\t400\t.\t+\t1\ttranscript_id "g1.t1"; gene_id "g1";\n'
)
r = m.augustus_gff_to_genbank(str(w / "genes.gtf"), str(w / "genome.fa"), 50, str(w / "train.gb"), verbose=False)
print(r["loci"], r["preview"][0].split()[1], r["preview"][3].strip())
```

- **Example output:**

```text
1 chr1_51-450 CDS             join(51..150,251..350)
```

## 17) `dnakit.annotation.augustus_gtf_to_gff()`

- **Purpose:** Convert AUGUSTUS GTF gene models to GFF (``gtf2gff.pl``); ``gff3=True`` writes GFF3 with ``ID``/``Parent`` attributes.
- **Method:** Runs `ID`, `Parent`, `gtf2gff.pl` through the DNAKit bounded executor.
- **API:** `dnakit.annotation.augustus_gtf_to_gff(gtf_path, output_path, *, gff3=False, print_exons=False, print_utrs=False, print_introns=False, include_stop_in_cds=False, augustus_bin_dir=None, verbose=True)`.
- **Input:** Required: `gtf_path`, `output_path`. Optional keyword-only: `gff3`, `print_exons`, `print_utrs`, `print_introns`, `include_stop_in_cds`, `augustus_bin_dir` (including `augustus_bin_dir`, `config_path` where applicable). Input files must exist; output paths must not.
- **Output:** A dictionary with `status`, `backend`, `version`, `command`, `return_code`, `outputs`, `log_tail`, `provenance` and parsed counts such as `lines`, `records`, `genes` or `transcripts`.
- **Example code:**

```python
import random, tempfile
from pathlib import Path
from dnakit.annotation import augustus_p04_conversion_functions as m
w = Path(tempfile.mkdtemp())
(w / "genes.gtf").write_text(
    'chr1\tAUGUSTUS\tgene\t101\t400\t.\t+\t.\tg1\n'
    'chr1\tAUGUSTUS\ttranscript\t101\t400\t.\t+\t.\tg1.t1\n'
    'chr1\tAUGUSTUS\tstart_codon\t101\t103\t.\t+\t0\ttranscript_id "g1.t1"; gene_id "g1";\n'
    'chr1\tAUGUSTUS\tCDS\t101\t200\t.\t+\t0\ttranscript_id "g1.t1"; gene_id "g1";\n'
    'chr1\tAUGUSTUS\tCDS\t301\t400\t.\t+\t1\ttranscript_id "g1.t1"; gene_id "g1";\n'
)
r = m.augustus_gtf_to_gff(str(w / "genes.gtf"), str(w / "genes.gff3"), gff3=True, print_exons=True, verbose=False)
print(r["preview"][0])
print(r["preview"][1])
```

- **Example output:**

```text
chr1	AUGUSTUS	gene	101	400	.	+	.	ID=g1;
chr1	AUGUSTUS	mRNA	101	400	.	+	.	ID=g1.t1;Parent=g1;
```

## 18) `dnakit.annotation.augustus_gth_to_gtf()`

- **Purpose:** Convert GenomeThreader GFF3 alignments to GTF (``gth2gtf.pl``).
- **Method:** Runs `gth2gtf.pl` through the DNAKit bounded executor.
- **API:** `dnakit.annotation.augustus_gth_to_gtf(gth_gff3_path, output_path, *, augustus_bin_dir=None, verbose=True)`.
- **Input:** Required: `gth_gff3_path`, `output_path`. Optional keyword-only: `augustus_bin_dir` (including `augustus_bin_dir`, `config_path` where applicable). Input files must exist; output paths must not.
- **Output:** A dictionary with `status`, `backend`, `version`, `command`, `return_code`, `outputs`, `log_tail`, `provenance` and parsed counts such as `lines`, `records`, `genes` or `transcripts`.
- **Example code:** not shown; the example needs data that is too large or an external tool that is not part of the AUGUSTUS package (see the tests and the manifest notes).

## 19) `dnakit.annotation.augustus_mask_n_regions()`

- **Purpose:** Write a GFF-like list of the N (gap) regions of a genome FASTA (``maskNregions.pl``).
- **Method:** Runs `maskNregions.pl` through the DNAKit bounded executor.
- **API:** `dnakit.annotation.augustus_mask_n_regions(genome_path, output_path, *, augustus_bin_dir=None, verbose=True)`.
- **Input:** Required: `genome_path`, `output_path`. Optional keyword-only: `augustus_bin_dir` (including `augustus_bin_dir`, `config_path` where applicable). Input files must exist; output paths must not.
- **Output:** A dictionary with `status`, `backend`, `version`, `command`, `return_code`, `outputs`, `log_tail`, `provenance` and parsed counts such as `lines`, `records`, `genes` or `transcripts`.
- **Example code:**

```python
import random, tempfile
from pathlib import Path
from dnakit.annotation import augustus_p04_conversion_functions as m
w = Path(tempfile.mkdtemp())
(w / "n.fa").write_text(">s1\nACGTNNNNACGT\n")
r = m.augustus_mask_n_regions(str(w / "n.fa"), str(w / "n.gff"), verbose=False)
print(r["preview"])
```

- **Example output:**

```text
['s1\tmaskN\tN-region\t5\t8\t.\t.\t.\tN-region']
```

## 20) `dnakit.annotation.augustus_merge_masking()`

- **Purpose:** Soft-mask every base that is masked in either of two versions of the same assembly (``merge_masking.pl``).
- **Method:** Runs `merge_masking.pl` through the DNAKit bounded executor.
- **API:** `dnakit.annotation.augustus_merge_masking(softmasked_path, second_softmasked_path, output_path, *, augustus_bin_dir=None, verbose=True)`.
- **Input:** Required: `softmasked_path`, `second_softmasked_path`, `output_path`. Optional keyword-only: `augustus_bin_dir` (including `augustus_bin_dir`, `config_path` where applicable). Input files must exist; output paths must not.
- **Output:** A dictionary with `status`, `backend`, `version`, `command`, `return_code`, `outputs`, `log_tail`, `provenance` and parsed counts such as `lines`, `records`, `genes` or `transcripts`.
- **Example code:**

```python
import random, tempfile
from pathlib import Path
from dnakit.annotation import augustus_p04_conversion_functions as m
w = Path(tempfile.mkdtemp())
rng = random.Random(7)
(w / "genome.fa").write_text(">chr1\n" + "".join(rng.choices("ACGT", k=4000)) + "\n")
r = m.augustus_merge_masking(str(w / "genome.fa"), str(w / "genome.fa"), str(w / "merged.fa"), verbose=False)
print(r["records"])
```

- **Example output:**

```text
1
```

## 21) `dnakit.annotation.augustus_partition_gene_models()`

- **Purpose:** Partition single-exon CDS gene models into gene clusters and write them as a GenBank file, one cluster per LOCUS (``partition_gtf2gb.pl``).
- **Method:** Runs `partition_gtf2gb.pl` through the DNAKit bounded executor.
- **API:** `dnakit.annotation.augustus_partition_gene_models(genome_path, gtf_path, output_path, *, min_cluster_size=None, max_flanking_region_size=None, augustus_bin_dir=None, verbose=True)`.
- **Input:** Required: `genome_path`, `gtf_path`, `output_path`. Optional keyword-only: `min_cluster_size`, `max_flanking_region_size`, `augustus_bin_dir` (including `augustus_bin_dir`, `config_path` where applicable). Input files must exist; output paths must not.
- **Output:** A dictionary with `status`, `backend`, `version`, `command`, `return_code`, `outputs`, `log_tail`, `provenance` and parsed counts such as `lines`, `records`, `genes` or `transcripts`.
- **Example code:** not shown; the example needs data that is too large or an external tool that is not part of the AUGUSTUS package (see the tests and the manifest notes).

## 22) `dnakit.annotation.augustus_random_split()`

- **Purpose:** Randomly split a GenBank gene set into a test set of ``test_size`` loci and a training set (``randomSplit.pl``).
- **Method:** Runs `randomSplit.pl`, `test_size` through the DNAKit bounded executor.
- **API:** `dnakit.annotation.augustus_random_split(genbank_path, test_size, train_path, test_path, *, augustus_bin_dir=None, verbose=True)`.
- **Input:** Required: `genbank_path`, `test_size`, `train_path`, `test_path`. Optional keyword-only: `augustus_bin_dir` (including `augustus_bin_dir`, `config_path` where applicable). Input files must exist; output paths must not.
- **Output:** A dictionary with `status`, `backend`, `version`, `command`, `return_code`, `outputs`, `log_tail`, `provenance` and parsed counts such as `lines`, `records`, `genes` or `transcripts`.
- **Example code:**

```python
import random, tempfile
from pathlib import Path
from dnakit.annotation import augustus_p04_conversion_functions as m
w = Path(tempfile.mkdtemp())
rng = random.Random(7)
(w / "genome.fa").write_text(">chr1\n" + "".join(rng.choices("ACGT", k=4000)) + "\n")
(w / "genes.gtf").write_text(
    'chr1\tAUGUSTUS\tgene\t101\t400\t.\t+\t.\tg1\n'
    'chr1\tAUGUSTUS\ttranscript\t101\t400\t.\t+\t.\tg1.t1\n'
    'chr1\tAUGUSTUS\tstart_codon\t101\t103\t.\t+\t0\ttranscript_id "g1.t1"; gene_id "g1";\n'
    'chr1\tAUGUSTUS\tCDS\t101\t200\t.\t+\t0\ttranscript_id "g1.t1"; gene_id "g1";\n'
    'chr1\tAUGUSTUS\tCDS\t301\t400\t.\t+\t1\ttranscript_id "g1.t1"; gene_id "g1";\n'
)
gb = m.augustus_gff_to_genbank(str(w / "genes.gtf"), str(w / "genome.fa"), 50, str(w / "train.gb"), verbose=False)
r = m.augustus_random_split(gb["outputs"]["output"], 1, str(w / "train.split.gb"), str(w / "test.split.gb"), verbose=False)
print(r["train_loci"] + r["test_loci"])
```

- **Example output:**

```text
1
```

## 23) `dnakit.annotation.augustus_scipio_gff_to_gff()`

- **Purpose:** Convert a Scipio GFF to an AUGUSTUS GFF that ``augustus_gff_to_genbank`` can read (``scipiogff2gff.pl``).
- **Method:** Runs `augustus_gff_to_genbank`, `scipiogff2gff.pl` through the DNAKit bounded executor.
- **API:** `dnakit.annotation.augustus_scipio_gff_to_gff(scipio_gff_path, output_path, *, augustus_bin_dir=None, verbose=True)`.
- **Input:** Required: `scipio_gff_path`, `output_path`. Optional keyword-only: `augustus_bin_dir` (including `augustus_bin_dir`, `config_path` where applicable). Input files must exist; output paths must not.
- **Output:** A dictionary with `status`, `backend`, `version`, `command`, `return_code`, `outputs`, `log_tail`, `provenance` and parsed counts such as `lines`, `records`, `genes` or `transcripts`.
- **Example code:** not shown; the example needs data that is too large or an external tool that is not part of the AUGUSTUS package (see the tests and the manifest notes).

## 24) `dnakit.annotation.augustus_scipio_yaml_to_gff()`

- **Purpose:** Convert Scipio YAML alignments to GFF (``yaml2gff.1.4.pl``); ``filter_status`` drops alignments with that status (for example ``incomplete``).
- **Method:** Runs `filter_status`, `incomplete`, `yaml2gff.1.4.pl` through the DNAKit bounded executor.
- **API:** `dnakit.annotation.augustus_scipio_yaml_to_gff(scipio_yaml_path, output_path, *, filter_status=None, augustus_bin_dir=None, verbose=True)`.
- **Input:** Required: `scipio_yaml_path`, `output_path`. Optional keyword-only: `filter_status`, `augustus_bin_dir` (including `augustus_bin_dir`, `config_path` where applicable). Input files must exist; output paths must not.
- **Output:** A dictionary with `status`, `backend`, `version`, `command`, `return_code`, `outputs`, `log_tail`, `provenance` and parsed counts such as `lines`, `records`, `genes` or `transcripts`.
- **Example code:** not shown; the example needs data that is too large or an external tool that is not part of the AUGUSTUS package (see the tests and the manifest notes).

## 25) `dnakit.annotation.augustus_summarize_base_content()`

- **Purpose:** Summarise the A, C, G, T, N and other base counts of each sequence (``summarizeACGTcontent.pl``).
- **Method:** Runs `summarizeACGTcontent.pl` through the DNAKit bounded executor.
- **API:** `dnakit.annotation.augustus_summarize_base_content(fasta_path, output_path, *, augustus_bin_dir=None, verbose=True)`.
- **Input:** Required: `fasta_path`, `output_path`. Optional keyword-only: `augustus_bin_dir` (including `augustus_bin_dir`, `config_path` where applicable). Input files must exist; output paths must not.
- **Output:** A dictionary with `status`, `backend`, `version`, `command`, `return_code`, `outputs`, `log_tail`, `provenance` and parsed counts such as `lines`, `records`, `genes` or `transcripts`.
- **Example code:**

```python
import random, tempfile
from pathlib import Path
from dnakit.annotation import augustus_p04_conversion_functions as m
w = Path(tempfile.mkdtemp())
rng = random.Random(7)
(w / "genome.fa").write_text(">chr1\n" + "".join(rng.choices("ACGT", k=4000)) + "\n")
r = m.augustus_summarize_base_content(str(w / "genome.fa"), str(w / "acgt.txt"), verbose=False)
print(r["preview"][-2])
```

- **Example output:**

```text
total 4000bp
```

