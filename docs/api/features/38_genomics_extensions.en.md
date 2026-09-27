# Genomics and DNA engineering

APIs for genotype PCA, ANI, breeding prediction, prokaryotic gene finding, methylation summaries, DNA storage and nanostructure export. See the [parameter FAQ](../../faq/api/38_genomics_extensions.md) for dependencies, full signatures and limitations, and the [source FAQ](../../faq/genomics_sources.md) for authors, papers and licenses.

<span id="pip-genomics"></span>
Install NumPy, Pyrodigal, pysam and PyFastANI together from this checkout with `python -m pip install -e ".[genomics]"`. Genotype PCA, ANI, genomic relationships, breeding prediction and prokaryotic CDS prediction below, plus [pysam-based BAM statistics](32_sequencing_analysis.md#bio-14), run inside Python without separately installed scientific command-line programs.

This extra does not include DNA storage, nanostructure dependencies or every external analysis. See the [43-analysis migration FAQ](../../faq/external_analysis_migration.md). These additions have not yet been published to PyPI.

## 1) Genomic relationships and breeding prediction

- **Purpose:** Calculate additive genomic relationships between samples and predict new samples' phenotypes from training genotypes and phenotypes.
- **Method:** VanRaden method 1 for relationships; fixed-lambda RR-BLUP/ridge regression with an intercept for prediction. Prediction uses training statistics alone for imputation, marker filtering and centering.
- **API:** `dnakit.genomics.genomic_relationship_matrix(genotypes[required], missing[optional])`; `dnakit.genomics.genomic_prediction(train_genotypes[required], phenotypes[required], predict_genotypes[required], ridge_lambda[required], missing[optional])`
- **Input:** Diploid dosage matrices with sample rows, marker columns and values in 0–2. Training and prediction must share marker order and allele orientation. Supply one finite phenotype per training sample and a positive penalty.
- **Output:** Returns relationship-matrix and genomic-prediction mappings with sample counts and predictions.
- **Example code:**

```python
from dnakit.genomics import genomic_relationship_matrix, genomic_prediction

genotypes = [[0, 0, 2], [1, 2, 1], [2, 1, 0], [0, 1, 1]]
relationship = genomic_relationship_matrix(genotypes)
prediction = genomic_prediction(
    genotypes,
    [2.0, 5.0, 4.0, 3.0],
    [[1, 1, 1]],
    ridge_lambda=1.0,
)
print(relationship["n_samples"])
print([round(value, 3) for value in prediction["predictions"]])
```

- **Example result:**

```text
4
[3.592]
```

[Parameters, missing values and validation](../../faq/api/38_genomics_extensions.md#breeding)

## 2) Prokaryotic gene prediction

- **Purpose:** Identify CDSs in a bacterial or archaeal linear contig, returning coordinates, coding sequences, proteins and model scores.
- **Method:** Call local Pyrodigal/Prodigal with upstream pretrained models by default. DNAKit converts coordinates to 0-based half-open intervals and records the backend version.
- **API:** `dnakit.genomics.predict_prokaryotic_genes(sequence[required], sequence_id[optional], meta[optional], training_sequence[optional])`
- **Input:** One nonempty, ungapped IUPAC DNA sequence. Ambiguous bases become N without shifting positions. The artificial sequence below demonstrates the interface and does not guarantee a gene call.
- **Output:** Returns a `GenePredictionResult`; the example shows `sequence_id`, `sequence_length`, `coordinate_system`.
- **Example code:**

```python
from dnakit.genomics import predict_prokaryotic_genes

sequence = "ATG" + "GCT" * 120 + "TAA"
result = predict_prokaryotic_genes(sequence, sequence_id="example", meta=True)
print(result.sequence_id, result.sequence_length)
print(result.coordinate_system)
```

- **Example result:** Sequence identity, length and coordinate convention from the result; gene predictions are stored in `result.genes`.

```text
example 366
0-based half-open
```

[Parameters, training mode and output fields](../../faq/api/38_genomics_extensions.md#genes)

## 3) Methylation counts

- **Purpose:** Read bedMethyl site counts and summarize methylation by chromosome, strand, modification code and motif.
- **Method:** Parse modkit's 18-column extended bedMethyl, filter records, then divide total modified counts by total valid coverage within each group. Keep 5mC and 5hmC separate.
- **API:** `dnakit.genomics.read_bedmethyl(path[required])`; `dnakit.genomics.summarize_methylation(records[required], modification_codes[optional], min_coverage[optional])`
- **Input:** A bedMethyl path, including gzip/BGZF files; column 11 is a percentage in 0–100. This example creates one temporary record with coverage 10 and modified count 5.
- **Output:** `read_bedmethyl` returns a record list; `summarize_methylation` returns modification fractions and grouped statistics.
- **Example code:**

```python
from pathlib import Path
from tempfile import TemporaryDirectory
from dnakit.genomics import read_bedmethyl, summarize_methylation

with TemporaryDirectory() as directory:
    path = Path(directory) / "example.bedmethyl"
    path.write_text(
        "chr1\t0\t1\tm\t10\t+\t0\t1\t255,0,0\t10\t50\t5\t5\t0\t0\t0\t0\t0\n",
        encoding="utf-8",
    )
    records = read_bedmethyl(path)
summary = summarize_methylation(records, modification_codes=["m"], min_coverage=10)
group = summary["groups"][0]
print(summary["records_selected"])
print(group["chromosome"], group["modification_code"], group["fraction_modified"])
```

- **Example result:**

```text
1
chr1 m 0.5
```

[Parameters, reader limits and statistical scope](../../faq/api/38_genomics_extensions.md#methylation)

## 4) DNA storage

- **Purpose:** Encode bytes as DNA segments and recover the original bytes from an intact packet.
- **Method:** Use Chamaeleo's fixed-Huffman Goldman codec, reorder segments by index during decoding, and check lengths and the SHA-256 digest.
- **API:** `dnakit.genomics.encode_dna_storage(data_hex[required], segment_bytes[optional])`; `dnakit.genomics.decode_dna_storage(packet[required])`
- **Input:** A hex-encoded byte string for encoding; the complete JSON packet, including DNA segments, indices, lengths and digest, for decoding.
- **Output:** Encoding returns a DNA storage package; decoding returns restored `data_hex`.
- **Example code:**

```python
from dnakit.genomics import encode_dna_storage, decode_dna_storage

packet = encode_dna_storage("DNAKit".encode().hex(), segment_bytes=32)
restored = decode_dna_storage(packet)
print(bytes.fromhex(restored["data_hex"]).decode())
```

- **Example result:**

```text
DNAKit
```

[Parameters, packet preservation and integrity checks](../../faq/api/38_genomics_extensions.md#storage)

## 5) DNA nanostructures and oxDNA export

- **Purpose:** Construct a complementary DNA duplex in scadnano and export oxDNA configuration and topology text.
- **Method:** Scadnano creates antiparallel complementary strands, assigns bases and converts coordinates. DNAKit validates sequences and geometry before export; it does not run folding optimization or molecular dynamics.
- **API:** `dnakit.genomics.design_dna_duplex(sequence[required])`; `dnakit.genomics.export_oxdna_design(design[required])`
- **Input:** A complete A/C/G/T sequence or scadnano JSON with complete sequence assignments on every strand. Export currently supports linear strands only.
- **Output:** Design returns a duplex-structure mapping; oxDNA export returns topology and coordinate text.
- **Example code:**

```python
from dnakit.genomics import design_dna_duplex, export_oxdna_design

design = design_dna_duplex("AACCGGTT")
export = export_oxdna_design(design)
print(export["strand_count"], export["nucleotide_count"])
print(export["topology"].splitlines()[0])
```

- **Example result:** Strand and nucleotide counts, followed by the oxDNA topology header's nucleotide and strand counts.

```text
2 16
16 2
```

[Parameters, saving files and export limitations](../../faq/api/38_genomics_extensions.md#nanotechnology)

## 6) Genotype PCA {#pca}

- **Purpose:** Calculate sample components and marker loadings from a prepared diploid dosage matrix.
- **Method:** Impute marker means, remove all-missing/invariant markers, frequency-standardize dosages and compute an SVD.
- **API:** `dnakit.genomics.genotype_pca(genotypes, components=2, sample_ids=None, verbose=True)`.
- **Input:** Samples × markers with dosages in [0, 2] and `None` for missing values. Install NumPy through `.[genomics]` or `.[bio-analysis]`.
- **Output:** Returns a mapping; the example shows `eigenvalues`, `marker_indices`.
- **Example code:**

```python
from dnakit.genomics import genotype_pca

result = genotype_pca(
    [[0, 0], [1, 2], [2, 1]],
    components=2,
    sample_ids=["sample_A", "sample_B", "sample_C"],
)
print([round(value, 6) for value in result["eigenvalues"]])
print(result["marker_indices"])
```

- **Example result:**

```text
[3.0, 1.0]
[0, 1]
```

Results also contain scores, loadings, variance ratios and preprocessing metadata. This dense-matrix PCA does not read BED/VCF, prune LD or reproduce a complete PLINK pipeline. [Parameters and eigenvalue scaling](../../faq/api/38_genomics_extensions.md#pca).

## 7) Genome ANI {#ani}

- **Purpose:** Estimate average nucleotide identity in homologous fragments of two related prokaryotic genomes.
- **Method:** Call FastANI through PyFastANI inside Python, preserving contig boundaries and computing query→reference only.
- **API:** `dnakit.genomics.average_nucleotide_identity(query, reference, threads=1, verbose=True)`.
- **Input:** Two assembly FASTA paths (gzip supported), or two DNA strings. Install `.[genomics]`; no FastANI executable is required.
- **Output:** Returns a mapping; the example shows `status`, `ani_percent`, `matching_fragments`, `query_fragments`.
- **Example code:** Run after preparing your own assembly files:

```python
from dnakit.genomics import average_nucleotide_identity

result = average_nucleotide_identity("query.fna", "reference.fna", threads=4)
print(result["status"], result["ani_percent"])
print(result["matching_fragments"], result["query_fragments"])
```

- **Example result:** Reportable matches have `status="matched"` and percentage ANI. Missing estimates have `status="no_match"` and `ani_percent=None`. Queries shorter than a complete fragment cannot yield an estimate. Values depend on the actual inputs; identical synthetic genomes validate near 100%. [Input conventions, directionality and full parameters](../../faq/api/38_genomics_extensions.md#ani).
