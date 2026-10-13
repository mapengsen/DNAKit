# 96 G18j Review: pyensembl, gget, biothings_client and genomepy

This page records the G18j review of four packages. The review handles only the 88 ledger rows whose status was `other`: it wraps 5 functions and gives each of the other 83 rows a status with a one-line reason. The wrappers live in new `dnakit.databases` modules and import directly, for example `from dnakit.databases import ensembl_species_record`.

- Offline operation: none of the 5 wrappers touches the network. The URL functions only build or match strings; the species function reads pyensembl's built-in species registry; the genomepy function reads only the `README.txt` and `assembly_report.txt` saved in a genome folder, and checks that both exist before it calls genomepy, so it never triggers a download.
- gget and biothings_client have no offline DNA-level function in this group. `gget_review_functions.py` and `biothings_client_review_functions.py` export an empty `__all__`, and their reasons are recorded in the ledger.
- Input files in the examples are declared with the `Example file` label; the other examples run directly.

## Review results

| Package | other/cap before | New wrappers | Other statuses |
| --- | ---: | ---: | --- |
| pyensembl | 56 | 4 | family-variant 15, needs-data 4, not-dna 3, internal 30 |
| gget | 4 | 0 | not-dna 2, needs-network 1, cli 1 |
| biothings_client | 5 | 0 | internal 3, not-dna 2 |
| genomepy | 23 | 1 | family-variant 1, internal 18, not-dna 2, cli 1 |
| Total | 88 | 5 | family-variant 16, internal 51, not-dna 9, needs-data 4, needs-network 1, cli 2 |

A `family-variant` row is a module alias of an upstream function that a wrapper already covers; the ledger's `dnakit_name` column names that wrapper.

## 1) `ensembl_dated_release_urls`

- **Purpose:** Give the GTF, cDNA, peptide and reference DNA download addresses of one dated annotation release (the new Ensembl platform, such as `2026_04`).
- **Method:** Check the GCA or GCF accession, the provider name and the annotation date (YYYY_MM), then call pyensembl's `make_dated_release_urls` to build the addresses. `include_alt` changes only the GTF file name; `genome_fasta_mask` sets the reference DNA file name (none, soft or hard). No network request is made.
- **API:** `ensembl_dated_release_urls(assembly_accession[required], provider[required], annotation_date[required], *, include_alt[optional], genome_fasta_mask[optional], verbose[optional])`
- **Input:** An accession such as `GCA_000001405.29`; a provider such as `ensembl`; an annotation date such as `2026_04`; `include_alt` as a boolean, default False; `genome_fasta_mask` as none, soft or hard, default none.
- **Output:** A dictionary with the five addresses `directory`, `gtf`, `cdna`, `pep` and `genome_fasta`, the echoed inputs and `provenance`.

**Example code:**

```python
from dnakit.databases import ensembl_dated_release_urls

result = ensembl_dated_release_urls(
    "GCA_000001405.29", "ensembl", "2026_04", genome_fasta_mask="hard", verbose=False
)
print(result["directory"])
print(result["gtf"])
print(result["genome_fasta"])
```

**Example output:**

```text
https://ftp.ebi.ac.uk/pub/ensemblorganisms/GCA/000/001/405/29/ensembl/2026_04
https://ftp.ebi.ac.uk/pub/ensemblorganisms/GCA/000/001/405/29/ensembl/2026_04/geneset/genes.gtf.gz
https://ftp.ebi.ac.uk/pub/ensemblorganisms/GCA/000/001/405/29/ensembl/2026_04/genome/hardmasked.fa.bgz
```

## 2) `ensembl_dated_genome_fasta_mask`

- **Purpose:** Read the mask (none, soft or hard) from the reference DNA address of a dated release.
- **Method:** Match the end of the address against `unmasked.fa.bgz`, `softmasked.fa.bgz` or `hardmasked.fa.bgz`. When nothing matches, `genome_fasta_mask` is None and `matched` is False.
- **API:** `ensembl_dated_genome_fasta_mask(url[required], *, verbose[optional])`
- **Input:** An address string of at most 2048 characters.
- **Output:** `url`, `genome_fasta_mask` and `matched`.

**Example code:**

```python
from dnakit.databases import ensembl_dated_genome_fasta_mask

base = "https://ftp.ebi.ac.uk/pub/ensemblorganisms/GCA/000/001/405/29/ensembl/2026_04/genome"
print(ensembl_dated_genome_fasta_mask(f"{base}/hardmasked.fa.bgz", verbose=False)["genome_fasta_mask"])
print(ensembl_dated_genome_fasta_mask(f"{base}/cdna.fa.bgz", verbose=False)["matched"])
```

**Example output:**

```text
hard
False
```

## 3) `ensembl_genome_fasta_source_options`

- **Purpose:** Read the FASTA type (toplevel or primary_assembly) and the mask of an official Ensembl reference DNA address.
- **Method:** The address must be an https address on `ftp.ensembl.org` or `ftp.ensemblgenomes.ebi.ac.uk` and its file name must follow pyensembl's official pattern. `dna_sm` means soft, `dna_rm` means hard and `dna` means none. Any other address returns `canonical` as False with the type and mask set to None.
- **API:** `ensembl_genome_fasta_source_options(source[required], *, verbose[optional])`
- **Input:** An address string of at most 2048 characters.
- **Output:** `canonical`, `genome_fasta_type` and `genome_fasta_mask`.

**Example code:**

```python
from dnakit.databases import ensembl_genome_fasta_source_options

url = (
    "https://ftp.ensembl.org/pub/release-112/fasta/homo_sapiens/dna/"
    "Homo_sapiens.GRCh38.dna_sm.toplevel.fa.gz"
)
result = ensembl_genome_fasta_source_options(url, verbose=False)
print(result["canonical"], result["genome_fasta_type"], result["genome_fasta_mask"])
print(ensembl_genome_fasta_source_options("/data/genome.fa.gz", verbose=False)["canonical"])
```

**Example output:**

```text
True toplevel soft
False
```

## 4) `ensembl_species_record`

- **Purpose:** Return the latin name, synonyms, the Ensembl release range of each reference assembly, and the accession and provider of the dated annotation.
- **Method:** `by="name"` looks the species up by latin or common name (such as `human`) in pyensembl's species registry; `by="reference"` looks it up by reference assembly name (such as `GRCh38`, case-insensitive). An unknown name raises `ConfigurationError`.
- **API:** `ensembl_species_record(query[required], *, by[optional], verbose[optional])`
- **Input:** A species name or a reference assembly name; `by` is `name` or `reference`, default `name`.
- **Output:** `latin_name`, `synonyms`, `division`, `ensembl_genomes`, `reference_assemblies` (each item has `reference`, `first_release` and `last_release`) and `dated_releases`.

**Example code:**

```python
from dnakit.databases import ensembl_species_record

record = ensembl_species_record("human", verbose=False)
print(record["latin_name"], record["division"])
for item in record["reference_assemblies"]:
    print(item["reference"], item["first_release"], item["last_release"])
print(ensembl_species_record("GRCh38", by="reference", verbose=False)["latin_name"])
```

**Example output:**

```text
homo_sapiens vertebrates
GRCh37 55 75
GRCh38 76 116
NCBI36 54 54
homo_sapiens
```

## 5) `map_genome_contig_names`

- **Purpose:** Map the chromosome or contig names of an installed local genome to the UCSC, NCBI or Ensembl naming scheme.
- **Method:** Read the `README.txt` of the genome folder (for the source provider) and its `assembly_report.txt` (the report genomepy saves at install time). When either file is missing the call fails before genomepy runs, so nothing is downloaded. When source and target provider match, or the report has no names for the target, `available` is False. A target name of `na` is returned as None.
- **API:** `map_genome_contig_names(genome_name[required], target_provider[required], *, genomes_dir[required], verbose[optional])`
- **Input:** A genome folder name (such as `demo`), a target provider (`ucsc`, `ncbi` or `ensembl`), and the `genomes_dir` that holds the folder. The example files below go into a new `genomes/demo` directory.
- **Output:** `genome`, `source_provider`, `target_provider`, `available`, `reason`, `record_count`, and `records` sorted by source name (each item has `source_name` and `target_name`).

**Example file:** `demo_README.txt`

```text
name: demo
provider: NCBI
assembly_accession: GCF_000001405.40
tax_id: 9606
species: Homo sapiens
genome_file: demo.fa
```

**Example file:** `demo_assembly_report.txt`

```text
Sequence-Name	Sequence-Role	Assigned-Molecule	Assigned-Molecule-Location/Type	GenBank-Accn	Relationship	RefSeq-Accn	Assembly-Unit	Sequence-Length	UCSC-style-name
1	assembled-molecule	1	Chromosome	CM000663.2	=	NC_000001.11	Primary Assembly	248956422	chr1
2	assembled-molecule	2	Chromosome	CM000664.2	=	NC_000002.12	Primary Assembly	242193529	chr2
MT	assembled-molecule	MT	Mitochondrion	J01415.2	=	NC_012920.1	non-nuclear	16569	chrM
KI270706.1	unplaced-scaffold	na	Unplaced-Scaffold	KI270706.1	=	na	Unlocalized	175055	na
```

**Example code:**

```python
import os
import shutil

from dnakit.databases import map_genome_contig_names

os.makedirs("genomes/demo", exist_ok=True)
shutil.copy("demo_README.txt", "genomes/demo/README.txt")
shutil.copy("demo_assembly_report.txt", "genomes/demo/assembly_report.txt")
result = map_genome_contig_names("demo", "UCSC", genomes_dir="genomes", verbose=False)
print(result["source_provider"], result["target_provider"], result["record_count"])
for record in result["records"]:
    print(record["source_name"], record["target_name"])
```

**Example output:**

```text
ncbi ucsc 4
1 chr1
2 chr2
KI270706.1 None
MT chrM
```

## Functions left unwrapped

- `needs-data`, 4 rows: `cached_release` and `genome_for_reference_name`, each with two module aliases, which return EnsemblRelease objects that depend on a local annotation cache.
- `needs-network`, 1 row: gget's `cellxgene`, which queries CELLxGENE over the network and needs the `cellxgene_census` package (not installed).
- `not-dna`, 9 rows: two species-name normalizers, one parser of annotation-date listings, two provider lists, gget's cBioPortal search and OpenAI text generation, and two biothings_client generic joins that combine two network clients.
- `cli`, 2 rows: gget's `setup` (an installer) and genomepy's `generate_env` (a shell export writer).
- `internal`, 51 rows: cache path and cache status helpers, destructive cleanup (`prune_genome_fastas`, `clean`), release-name validators, genomepy README reading and writing, client factories, a request cache helper and configuration writing.
- `family-variant`, 16 rows: module aliases of the functions above; see each wrapper section.

## Tests

The test file is `others/tests/unit/test_review_G18j.py`, with 12 cases, all offline. Every wrapper has at least one case that compares its output with a hand-built string, a direct upstream call or a known registry value. The missing-report case of genomepy replaces the upstream download function and asserts that it is never called.
