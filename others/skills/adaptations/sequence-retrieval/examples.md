# Sequence Retrieval Examples

Each Python block is independent and uses DNAKit's public APIs. Run in a working
directory where output filenames do not already exist. Network access is required.
FASTA/NCBI GenBank/ENA EMBL need only the base package. Format conversion in
Example 5 requires `pip install 'dnakit[bio-formats]'`.

## Example 1: Find E. coli K-12 Genome

```python
from pathlib import Path
import dnakit

result = dnakit.search.ncbi.entrez(
    database="nuccore",
    organism="Escherichia coli",
    strain="K-12",
    seq_type="complete_genome",
    retmax=3,
)
accessions = dnakit.search.ncbi.accessions(uids=result.metadata["uids"])
print("Search accessions:", accessions)

# An explicit reference choice, independent of the returned search ranking.
download = dnakit.download.ncbi.sequence(
    accession="NC_000913.3",
    format="fasta",
    output_path="ecoli_k12.fasta",
)
fasta = Path(download.files[0].path).read_text(encoding="utf-8")
length = sum(len(line.strip()) for line in fasta.splitlines() if not line.startswith(">"))
print(f"Genome size: {length} bases")
print("Source and hashes:", download.manifest_path)
```

## Example 2: Get Human BRCA1 Transcripts

```python
import dnakit

result = dnakit.search.ncbi.entrez(
    "nuccore", organism="Homo sapiens", gene="BRCA1",
    seq_type="mrna", source="refseq", retmax=5,
)
print(f"Found {result.total_count} matching records; showing {len(result.records)}")
accessions = dnakit.search.ncbi.accessions(result.metadata["uids"])
for record in result.records:
    print(record)
if not accessions:
    raise RuntimeError("No BRCA1 transcript candidates; review the query.")

# Demonstrates downloading one candidate. First search hit is not necessarily
# MANE Select or the canonical transcript; review before biological analysis.
candidate = accessions[0]
genbank = dnakit.download.ncbi.sequence(
    candidate, "brca1_candidate.gb", format="genbank",
)
print(candidate, genbank.files[0].path)
```

## Example 3: SARS-CoV-2 Reference Genome

```python
from pathlib import Path
import dnakit

# Verify the explicit version; do not assume it is the current/latest version.
reference = dnakit.search.ncbi.accession("NC_045512.2")
if not reference.records:
    raise RuntimeError("The requested reference was not found.")
print(reference.records)
genome = dnakit.download.ncbi.sequence("NC_045512.2", "sars_cov_2.fasta")
print(Path(genome.files[0].path).read_text(encoding="utf-8")[:200])
```

## Example 4: Compare RefSeq and Submitted Records

```python
import dnakit

# Search each source separately, so the top five hits need not contain both.
for source in ("refseq", "genbank"):
    result = dnakit.search.ncbi.entrez(
        "nuccore", organism="Escherichia coli", strain="K-12",
        seq_type="complete_genome", source=source, retmax=5,
    )
    accessions = dnakit.search.ncbi.accessions(result.metadata["uids"])
    print(source, accessions)
    print(result.records)

# NZ_ is also RefSeq. Non-RefSeq status does not prove ENA availability or
# experimental validation. Verify cross-references and exact strain/version.
```

## Example 5: Multi-Format Retrieval

```python
import dnakit

accession = "NC_000913.3"
fasta = dnakit.download.ncbi.sequence(accession, "ecoli_multiformat.fasta", format="fasta")
genbank = dnakit.download.ncbi.sequence(accession, "ecoli_multiformat.gb", format="genbank")

# Requires dnakit[bio-formats]. NCBI supplies GenBank; DNAKit uses Biopython
# to convert it to EMBL and records that conversion in the manifest.
embl = dnakit.download.ncbi.sequence(accession, "ecoli_multiformat.embl", format="embl")
print(fasta.files[0].path, genbank.files[0].path, embl.files[0].path)
```

## Example 6: ENA Metadata and Native Formats

```python
import dnakit

# This is an INSDC accession, not the RefSeq NC_000913.3 identifier.
metadata = dnakit.search.ena.ena_search(
    result="sequence", query='accession="U00096"',
    fields=("accession", "description", "scientific_name", "tax_id"), limit=1,
)
print(metadata.records)
fasta = dnakit.download.ena.sequence("U00096.3", "ena_ecoli.fasta")
embl = dnakit.download.ena.sequence("U00096.3", "ena_ecoli.embl", format="embl")
print(fasta.files[0].path, embl.files[0].path)
```

## Example 7: Protein Sequence and Annotations

```python
import dnakit

result = dnakit.search.ncbi.entrez(
    "protein", organism="Homo sapiens", gene="BRCA1", source="refseq", retmax=3,
)
accessions = dnakit.search.ncbi.accessions(result.metadata["uids"], database="protein")
if not accessions:
    raise RuntimeError("No protein candidates; review the query.")
print(result.records)
candidate = accessions[0]
fasta = dnakit.download.ncbi.sequence(candidate, "brca1_protein.fasta", database="protein")
genpept = dnakit.download.ncbi.sequence(
    candidate, "brca1_protein.gp", database="protein", format="genbank",
)
print(fasta.files[0].path, genpept.files[0].path)
```
