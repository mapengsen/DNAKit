# Gene names and intervals from genomic DNA

`dnakit.annotation.annotate_genes()` maps genomic DNA to a reference and transfers named gene intervals from its GFF3. Supply a species name or NCBI TaxID; `assembly` is optional. On the first request without an assembly, DNAKit selects the species' current RefSeq reference assembly. The exact downloaded assembly and annotation are then pinned in the project cache.

This transfers existing reference annotations. It does not perform de novo gene prediction, name novel genes, or annotate spliced cDNA/transcripts. Missing gene names remain null; stable IDs or locus tags remain available.

## Command line

Install the working copy:

```bash
python -m pip install -e .
```

Install the external mapper:

```bash
conda install -c conda-forge -c bioconda minimap2
```

Annotate a single-record genomic FASTA using an automatically selected assembly:

```bash
python -m dnakit.cli annotate-genes \
  --input query.fasta \
  --species "Homo sapiens" \
  --cache-dir .dnakit/gene_annotation \
  --output genes.json
```

Specify an exact assembly:

```bash
python -m dnakit.cli annotate-genes \
  --input query.fasta \
  --species "Homo sapiens" \
  --assembly GCF_000001405.40 \
  --cache-dir .dnakit/gene_annotation \
  --output genes.json
```

Run offline with a matched local reference/annotation pair:

```bash
python -m dnakit.cli annotate-genes \
  --input query.fasta \
  --species "Homo sapiens" \
  --reference-fasta reference.fna \
  --annotation-gff annotation.gff3 \
  --cache-dir .dnakit/gene_annotation \
  --output genes.json \
  --no-download
```

Use `--minimap2 /path/to/minimap2` for an explicit executable, `--threads 4` for concurrency, and `--overwrite` to replace an existing output. FASTA/GFF3 may be gzip-compressed. Each invocation processes one query; Python callers can loop over records and reuse the cache.

## Python

```python
from dnakit.annotation import GeneAnnotationConfig, annotate_genes
from dnakit.io import read

with read("query.fasta", format="fasta") as source:
    record = next(iter(source))

result = annotate_genes(
    record.sequence.symbols,
    species="Homo sapiens",
    assembly=None,
    sequence_id=record.id,
    config=GeneAnnotationConfig(cache_dir=".dnakit/gene_annotation", threads=4),
)
print(result.reference["assembly"])
print(result.to_dict())
```

Assembly requests accept versioned GCA_/GCF_ accessions, existing fixed aliases, or unambiguous species-scoped assembly names. Explicit old accessions are not replaced with a newer version. Multiple matches (including paired GenBank/RefSeq assemblies sharing a name) require an exact accession. Missing references or gene annotations produce actionable errors.

## Results

- `reference`: actual assembly, annotation metadata, file checksums, selection method, and species-validation status.
- `gene_records`: overlapping gene instances, names/IDs, strand, reference intervals, and projected query parts. A missing `gene_name` is null; `display_name` falls back to locus tag/ID.
- `segments`: a disjoint partition of the whole query into `gene`, `intergenic`, and `unknown`. A segment can refer to several overlapping genes.
- `alignments`: mapping coordinates, strand, MAPQ, identity and acceptance status.

Python `Interval` objects are **0-based half-open**. JSON from `to_dict()` and the CLI is **1-based closed** throughout. A gene record's start/end is an envelope; `parts` describes the confidently mapped bases and may have holes. Reverse strands and indels are projected using CIGAR blocks. Circular origin-spanning genes are supported when the reference is explicitly marked circular in GFF3.

For example, genes A:1–12000, B:18001–45000 and C:43000–52000 produce three gene records and these segments:

```text
1–12000      gene         A
12001–18000  intergenic
18001–42999  gene         B
43000–45000  gene         B + C
45001–52000  gene         C
```

## Reproducibility and scope

Identical species/assembly requests reuse a pinned reference without querying NCBI. SHA256 protects the cached FASTA, GFF3 and minimap2 index; indexes are keyed by reference content, preset and mapper version. Use a new cache directory to refresh reference/annotation selection. `--no-download` permits only pinned references or local file pairs.

Downloaded assemblies are checked against the requested taxonomy (including strain-to-species ancestry). Offline local species are caller-declared (`species_validated=false`). Local assembly declarations, sequence IDs and bounds are checked when available; callers remain responsible for providing a matching pair. GFF3 needs gene-level features, not only transcripts/CDS. Contig identifiers must match exactly; chromosome aliases are not guessed.

The default preset is `asm5` for long accurate genomic DNA. Defaults accept primary mappings with MAPQ ≥20, identity ≥0.90 and ≥100 aligned bases; MAPQ=255 is unknown. These are configurable filters, not calibrated accuracy guarantees.

`intergenic` means reliably mapped bases outside genes in the selected annotation. Unmapped/low-confidence regions, insertions, overlapping trusted primary mappings and contigs missing annotation coverage information are `unknown`. Absence of annotation is not evidence of biological absence of a gene.

Downloading a large reference and building a single-part index can require substantial storage/RAM. Configuration bounds reference/query length, downloads, GFF records, index/alignment size, and external process duration. Imports perform no downloads or external execution.

Run coordinate, reference resolution, cache and CLI tests; installing minimap2 also enables the real-process test:

```bash
PYTHONPATH=src python -m pytest -q others/tests/unit/test_gene_annotation.py
```

Sources: [NCBI reference selection](https://www.ncbi.nlm.nih.gov/datasets/docs/v2/data-processing/policies-annotation/genome-processing/reference-selection/), [genome package](https://www.ncbi.nlm.nih.gov/datasets/docs/v2/reference-docs/data-packages/genome/), [GFF3 conventions](https://www.ncbi.nlm.nih.gov/datasets/docs/v2/reference-docs/file-formats/annotation-files/about-ncbi-gff3/), [minimap2 manual](https://lh3.github.io/minimap2/minimap2.html).

This version adds multi-species gene annotation with optional assembly selection, overlapping intervals and reproducible reference caches.
