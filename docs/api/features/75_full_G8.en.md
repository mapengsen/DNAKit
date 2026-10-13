# pyensembl, genomepy and NCBI Datasets genome functions

These functions wrap DNA sequence, genomic interval, annotation-file and gene-identifier operations from [pyensembl](https://github.com/openvax/pyensembl), [genomepy](https://github.com/vanheeringen-lab/genomepy) and [ncbi-datasets-pylib](https://github.com/ncbi/datasets). Most of them read only the local files or records that the caller passes and download nothing; none performs alignment or genotype analysis. The NCBI Datasets assembly queries need the network, and local tests use a fake client. [pybiomart](https://github.com/jrowen/pybiomart) exposes no module-level public functions, so this group does not include it. All dependencies are optional, installed with the extras `bio-pyensembl`, `bio-genomepy` and `bio-ncbi-datasets`.

## Example data

**Example file:** `toy.fa`

```text
>chr1 test
ACGT
GG
>chr2
TTTT
```

Small FASTA with two sequences.

**Example file:** `toy.gtf`

```text
chr1	src	gene	10	100	.	+	.	gene_id "G1";
chr2	src	gene	5	50	.	-	.	gene_id "G2";
```

Two-line GTF annotation (tab-separated).

**Example file:** `toy_n.fa`

```text
>chr1
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTNNNNACGTACGTACGTACGTACGT
>chr2 desc
GGGGCCCCAAAATTTT
>scaf_9
NNNNACGTAC
```

FASTA with N runs: four Ns in the middle of chr1 and four at the start of scaf_9.

## 1) `dnakit.databases.read_fasta_records()`

- **Purpose:** Read a FASTA or gzipped FASTA file and return sequence records sorted by identifier.
- **Method:** Calls `pyensembl.fasta.parse_fasta_dictionary`. Rejects files larger than `max_bytes` (default 64 MB) before reading. No sequence statistics are computed.
- **API:** `dnakit.databases.read_fasta_records(fasta_path, *, max_bytes=64000000, verbose=True)`.
- **Input:** Path to a FASTA file (`.gz` allowed); `max_bytes` is the file-size limit.
- **Output:** A dictionary with `records` (`id`, `length`, `sequence`), `record_count` and `provenance`.
- **Example code:**

```python
from dnakit.databases import read_fasta_records

result = read_fasta_records("toy.fa", verbose=False)
print(result["record_count"])
for record in result["records"]:
    print(record["id"], record["length"], record["sequence"])
```

- **Example output:**

```text
2
chr1 6 ACGTGG
chr2 4 TTTT
```

## 2) `dnakit.databases.read_first_fasta_id()`

- **Purpose:** Return the identifier of the first FASTA record, or `None` when the file has no records.
- **Method:** Calls `pyensembl.fasta.first_fasta_id`, reading only up to the first `>` header line.
- **API:** `dnakit.databases.read_first_fasta_id(fasta_path, *, verbose=True)`.
- **Input:** Path to a FASTA file (`.gz` allowed).
- **Output:** A dictionary with `identifier` (a string or `None`).
- **Example code:**

```python
from dnakit.databases import read_first_fasta_id

print(read_first_fasta_id("toy.fa", verbose=False)["identifier"])
```

- **Example output:**

```text
chr1
```

## 3) `dnakit.databases.merge_genomic_intervals()`

- **Purpose:** Merge overlapping or adjacent 1-based inclusive intervals.
- **Method:** Calls `pyensembl.common.merge_intervals`. Input is sorted first; `(a, b)` and `(b + 1, c)` count as adjacent and are merged.
- **API:** `dnakit.databases.merge_genomic_intervals(intervals, *, max_intervals=1000000, verbose=True)`.
- **Input:** A list of `(start, end)` pairs, 1-based, with `start` not after `end`; `max_intervals` limits the count.
- **Output:** `intervals` (dictionaries with `start` and `end`, ascending), `input_count` and `provenance`.
- **Example code:**

```python
from dnakit.databases import merge_genomic_intervals

result = merge_genomic_intervals([(5, 9), (1, 3), (4, 4), (12, 15)], verbose=False)
print(result["input_count"])
for item in result["intervals"]:
    print(item["start"], item["end"])
```

- **Example output:**

```text
4
1 9
12 15
```

## 4) `dnakit.databases.find_nearest_genomic_locus()`

- **Purpose:** Find the locus nearest to a query interval.
- **Method:** Builds `pyensembl.locus.Locus` objects and calls `pyensembl.search.find_nearest_locus`. On a tie the earlier locus in the input is kept.
- **API:** `dnakit.databases.find_nearest_genomic_locus(start, end, loci, *, max_loci=100000, verbose=True)`.
- **Input:** Query `start` and `end` (1-based inclusive); `loci` is a list of dictionaries with `contig`, `start`, `end` and `strand`.
- **Output:** `distance` (`None` when there are no loci), `nearest_index`, `nearest` and `provenance`.
- **Example code:**

```python
from dnakit.databases import find_nearest_genomic_locus

loci = [
    {"contig": "1", "start": 100, "end": 200, "strand": "+"},
    {"contig": "2", "start": 500, "end": 600, "strand": "-"},
]
result = find_nearest_genomic_locus(300, 310, loci, verbose=False)
print(result["distance"], result["nearest_index"], result["nearest"]["contig"])
```

- **Example output:**

```text
100 0 1
```

## 5) `dnakit.databases.stable_ensembl_id()`

- **Purpose:** Strip the Ensembl version suffix such as `.15` from an identifier.
- **Method:** Calls `pyensembl.gene_name_aliases.stable_ensembl_gene_id`. Identifiers in other formats are returned unchanged.
- **API:** `dnakit.databases.stable_ensembl_id(identifier, *, verbose=True)`.
- **Input:** An identifier string.
- **Output:** `identifier` and `stable_id`.
- **Example code:**

```python
from dnakit.databases import stable_ensembl_id

print(stable_ensembl_id("ENSG00000139618.15", verbose=False)["stable_id"])
print(stable_ensembl_id("AT1G01010.1", verbose=False)["stable_id"])
```

- **Example output:**

```text
ENSG00000139618
AT1G01010.1
```

## 6) `dnakit.databases.normalize_chromosome_name()`

- **Purpose:** Normalize a chromosome or sequence name, for example `chrx` to `chrX` and the integer `7` to `"7"`.
- **Method:** Calls `pyensembl.normalization.normalize_chromosome`. Empty names and `"0"` raise `ConfigurationError`.
- **API:** `dnakit.databases.normalize_chromosome_name(name, *, verbose=True)`.
- **Input:** A string or an integer (booleans are rejected).
- **Output:** `input` and `chromosome`.
- **Example code:**

```python
from dnakit.databases import normalize_chromosome_name

for name in ["chrx", "x", 7, "chrUn_gl000212"]:
    print(normalize_chromosome_name(name, verbose=False)["chromosome"])
```

- **Example output:**

```text
chrX
X
7
chrUn_gl000212
```

## 7) `dnakit.databases.normalize_strand_symbol()`

- **Purpose:** Write a strand as `+` or `-`.
- **Method:** Calls `pyensembl.normalization.normalize_strand`, which accepts `+`, `-`, `1`, `-1` as strings and `1`, `-1` as integers.
- **API:** `dnakit.databases.normalize_strand_symbol(strand, *, verbose=True)`.
- **Input:** A string or an integer (booleans are rejected).
- **Output:** `input` and `strand`.
- **Example code:**

```python
from dnakit.databases import normalize_strand_symbol

print(normalize_strand_symbol("1", verbose=False)["strand"])
print(normalize_strand_symbol(-1, verbose=False)["strand"])
```

- **Example output:**

```text
+
-
```

## 8) `dnakit.databases.normalize_gene_aliases()`

- **Purpose:** Normalize a mapping from gene names to gene identifiers.
- **Method:** Calls `pyensembl.gene_name_aliases.normalize_aliases`. A single string ID counts as one value; IDs under each name are de-duplicated and sorted.
- **API:** `dnakit.databases.normalize_gene_aliases(aliases, *, max_names=100000, verbose=True)`.
- **Input:** A mapping from gene name to an ID string or a list of IDs; `max_names` limits the number of names.
- **Output:** `records` (`name`, `gene_ids`, sorted by name), `name_count` and `provenance`.
- **Example code:**

```python
from dnakit.databases import normalize_gene_aliases

result = normalize_gene_aliases({"BRCA2": "ENSG1.2", "X": ["b", "a", "a"]}, verbose=False)
for record in result["records"]:
    print(record["name"], record["gene_ids"])
```

- **Example output:**

```text
BRCA2 ['ENSG1.2']
X ['a', 'b']
```

## 9) `dnakit.databases.match_identifier_version()`

- **Purpose:** Find which installed form of an Ensembl identifier matches a given identifier.
- **Method:** Compares only stable IDs (version removed) first, then calls `pyensembl.versioned_ids.match_version`. A stable ID that is not installed returns `None`, so another gene is never matched.
- **API:** `dnakit.databases.match_identifier_version(identifier, installed, *, verbose=True)`.
- **Input:** An identifier; `installed` maps each installed identifier to its recorded version or `None`.
- **Output:** `matched` (the installed form or `None`), `installed_count`, `compared_count` and `provenance`.
- **Example code:**

```python
from dnakit.databases import match_identifier_version

print(match_identifier_version("ENSG1.2", {"ENSG1": 2}, verbose=False)["matched"])
print(match_identifier_version("ENSG9", {"ENSG1": 2}, verbose=False)["matched"])
```

- **Example output:**

```text
ENSG1
None
```

## 10) `dnakit.databases.lookup_sequence_version_fallback()`

- **Purpose:** Look up a sequence by identifier, tolerating `.N` version suffixes.
- **Method:** Calls `pyensembl.sequence_data.lookup_sequence_with_version_fallback`. A plain dictionary has no version index that pyensembl builds for FASTA data, so a bare ID does not match a versioned key.
- **API:** `dnakit.databases.lookup_sequence_version_fallback(sequences, identifier, *, version=None, verbose=True)`.
- **Input:** A `{identifier: sequence}` dictionary, an identifier and an optional `version`.
- **Output:** `found`, `sequence`, `identifier`, `version` and `provenance`.
- **Example code:**

```python
from dnakit.databases import lookup_sequence_version_fallback

sequences = {"ENSP1.3": "MKV", "ENST2": "ACGT"}
print(lookup_sequence_version_fallback(sequences, "ENST2.5", verbose=False)["sequence"])
print(lookup_sequence_version_fallback(sequences, "ENSP1.3", verbose=False)["sequence"])
```

- **Example output:**

```text
ACGT
MKV
```

## 11) `dnakit.databases.normalize_genome_reference_name()`

- **Purpose:** Write a reference genome name in its standard form, for example `grch38` as `GRCh38`.
- **Method:** Calls `pyensembl.reference_name.normalize_reference_name`, ignoring case. Names missing from the pyensembl registry raise `ConfigurationError`.
- **API:** `dnakit.databases.normalize_genome_reference_name(name, *, verbose=True)`.
- **Input:** A reference genome name string.
- **Output:** `input` and `reference`.
- **Example code:**

```python
from dnakit.databases import normalize_genome_reference_name

print(normalize_genome_reference_name("grch38", verbose=False)["reference"])
```

- **Example output:**

```text
GRCh38
```

## 12) `dnakit.databases.ensembl_reference_name()`

- **Purpose:** Return the reference assembly that a species uses in an Ensembl release.
- **Method:** Calls `pyensembl.reference_name.which_reference`, which reads only the registry bundled with pyensembl. No network access.
- **API:** `dnakit.databases.ensembl_reference_name(species_name, ensembl_release, *, verbose=True)`.
- **Input:** A species name such as `human`; an Ensembl release number (integer) or an annotation date such as `2026_04`.
- **Output:** `species`, `ensembl_release`, `reference` and `provenance`.
- **Example code:**

```python
from dnakit.databases import ensembl_reference_name

print(ensembl_reference_name("human", 75, verbose=False)["reference"])
print(ensembl_reference_name("human", 110, verbose=False)["reference"])
```

- **Example output:**

```text
GRCh37
GRCh38
```

## 13) `dnakit.databases.ensembl_gtf_url()`

- **Purpose:** Build the download URL of the gene annotation GTF for a species and Ensembl release.
- **Method:** Calls `pyensembl.ensembl_url_templates.make_gtf_url`. It only formats the string and makes no request.
- **API:** `dnakit.databases.ensembl_gtf_url(ensembl_release, species, *, verbose=True)`.
- **Input:** An Ensembl release number (integer) and a species name.
- **Output:** `url`, `species` and `ensembl_release`.
- **Example code:**

```python
from dnakit.databases import ensembl_gtf_url

print(ensembl_gtf_url(110, "human", verbose=False)["url"])
```

- **Example output:**

```text
https://ftp.ensembl.org/pub/release-110/gtf/homo_sapiens/Homo_sapiens.GRCh38.110.chr_patch_hapl_scaff.gtf.gz
```

## 14) `dnakit.databases.ensembl_genome_fasta_url()`

- **Purpose:** Build the download URL of the reference DNA FASTA.
- **Method:** Calls `pyensembl.ensembl_url_templates.make_genome_fasta_url`. `fasta_type` is `toplevel` or `primary_assembly`; `mask` is `none`, `soft` or `hard`. Only the string is built.
- **API:** `dnakit.databases.ensembl_genome_fasta_url(ensembl_release, species, *, fasta_type="toplevel", mask="none", verbose=True)`.
- **Input:** Ensembl release number, species name, `fasta_type` and `mask`.
- **Output:** `url` with the chosen `fasta_type` and `mask`.
- **Example code:**

```python
from dnakit.databases import ensembl_genome_fasta_url

result = ensembl_genome_fasta_url(110, "human", fasta_type="primary_assembly", mask="soft", verbose=False)
print(result["url"])
```

- **Example output:**

```text
https://ftp.ensembl.org/pub/release-110/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa.gz
```

## 15) `dnakit.databases.ensembl_sequence_fasta_url()`

- **Purpose:** Build the download URL of the cDNA, ncRNA or protein FASTA.
- **Method:** Calls `pyensembl.ensembl_url_templates.make_fasta_url`. `sequence_type` accepts only `cdna`, `ncrna` and `pep`. Only the string is built.
- **API:** `dnakit.databases.ensembl_sequence_fasta_url(ensembl_release, species, sequence_type, *, verbose=True)`.
- **Input:** Ensembl release number, species name and `sequence_type`.
- **Output:** `url`, `sequence_type`, `species` and `ensembl_release`.
- **Example code:**

```python
from dnakit.databases import ensembl_sequence_fasta_url

print(ensembl_sequence_fasta_url(110, "mouse", "pep", verbose=False)["url"])
```

- **Example output:**

```text
https://ftp.ensembl.org/pub/release-110/fasta/mus_musculus/pep/Mus_musculus.GRCm39.pep.all.fa.gz
```

## 16) `dnakit.databases.read_annotation_records()`

- **Purpose:** Read a GTF (9 columns) or BED12 (12 columns) annotation file and return one record per line.
- **Method:** Calls `genomepy.annotation.utils.read_annot`, choosing the format from the column count; other counts raise `ConfigurationError`. The file size is limited by `max_bytes` (default 64 MB).
- **API:** `dnakit.databases.read_annotation_records(annotation_path, *, max_bytes=64000000, max_records=1000000, verbose=True)`.
- **Input:** Path to the annotation file; `max_bytes` and `max_records` are limits.
- **Output:** `file_format` (`gtf` or `bed12`), `columns`, `record_count`, `records` and `provenance`.
- **Example code:**

```python
from dnakit.databases import read_annotation_records

result = read_annotation_records("toy.gtf", verbose=False)
print(result["file_format"], result["record_count"])
for record in result["records"]:
    print(record["seqname"], record["start"], record["end"], record["strand"])
```

- **Example output:**

```text
gtf 2
chr1 10 100 +
chr2 5 50 -
```

## 17) `dnakit.databases.write_annotation_records()`

- **Purpose:** Write GTF or BED12 records to a new tab-separated file without a header line.
- **Method:** Checks that each record has exactly the columns of the chosen format, then calls `genomepy.annotation.utils.write_annot`. The output file must not exist; existing files are never overwritten.
- **API:** `dnakit.databases.write_annotation_records(records, output_path, *, file_format="gtf", max_records=1000000, verbose=True)`.
- **Input:** A list of records; `output_path`; `file_format` (`gtf` or `bed12`).
- **Output:** `output_path`, `file_format`, `record_count` and `provenance`.
- **Example code:**

```python
from pathlib import Path

from dnakit.databases import write_annotation_records

Path("written.gtf").unlink(missing_ok=True)
record = {
    "seqname": "chr1",
    "source": "s",
    "feature": "gene",
    "start": 1,
    "end": 5,
    "score": ".",
    "strand": "+",
    "frame": ".",
    "attribute": 'gene_id "X";',
}
result = write_annotation_records([record], "written.gtf", verbose=False)
print(result["record_count"])
print(Path("written.gtf").read_text(encoding="utf-8"), end="")
```

- **Example output:**

```text
1
chr1	s	gene	1	5	.	+	.	gene_id "X";
```

## 18) `dnakit.databases.filter_annotation_contigs()`

- **Purpose:** Keep the annotation records whose contig name matches a regular expression.
- **Method:** Calls `genomepy.annotation.filter_regex`, which applies `re.match` to the first column. With `invert_match` true, the non-matching records are kept.
- **API:** `dnakit.databases.filter_annotation_contigs(records, pattern, *, file_format="gtf", invert_match=False, max_records=1000000, verbose=True)`.
- **Input:** A list of annotation records; `pattern` (a regular expression); `file_format`; `invert_match`.
- **Output:** `records` (kept records in their original order), `kept_count`, `removed_count` and `provenance`.
- **Example code:**

```python
from dnakit.databases import filter_annotation_contigs, read_annotation_records

records = read_annotation_records("toy.gtf", verbose=False)["records"]
kept = filter_annotation_contigs(records, "chr2", verbose=False)
print(kept["kept_count"], kept["removed_count"], kept["records"][0]["seqname"])
```

- **Example output:**

```text
1 1 chr2
```

## 19) `dnakit.databases.write_gap_bed()`

- **Purpose:** Write the runs of N (gaps) in a FASTA file to a BED file.
- **Method:** Calls `genomepy.genome.generate_gap_bed`, which reads the file through pyfaidx. Coordinates are 0-based half-open. pyfaidx writes a `.fai` index beside the input FASTA.
- **API:** `dnakit.databases.write_gap_bed(fasta_path, output_path, *, verbose=True)`.
- **Input:** Path to a FASTA file; a new output path.
- **Output:** `output_path`, `gap_count` and `provenance`.
- **Example code:**

```python
from pathlib import Path

from dnakit.databases import write_gap_bed

Path("gaps.bed").unlink(missing_ok=True)
result = write_gap_bed("toy_n.fa", "gaps.bed", verbose=False)
print(result["gap_count"])
print(Path("gaps.bed").read_text(encoding="utf-8"), end="")
```

- **Example output:**

```text
2
chr1	40	44
scaf_9	0	4
```

## 20) `dnakit.databases.write_fasta_sizes()`

- **Purpose:** Write the name and length of every sequence in a FASTA file.
- **Method:** Calls `genomepy.genome.generate_fa_sizes`, which writes a two-column, tab-separated sizes file.
- **API:** `dnakit.databases.write_fasta_sizes(fasta_path, output_path, *, verbose=True)`.
- **Input:** Path to a FASTA file; a new output path.
- **Output:** `output_path`, `sequence_count` and `provenance`.
- **Example code:**

```python
from pathlib import Path

from dnakit.databases import write_fasta_sizes

Path("sizes.txt").unlink(missing_ok=True)
result = write_fasta_sizes("toy_n.fa", "sizes.txt", verbose=False)
print(result["sequence_count"])
print(Path("sizes.txt").read_text(encoding="utf-8"), end="")
```

- **Example output:**

```text
3
chr1	64
chr2	16
scaf_9	10
```

## 21) `dnakit.databases.filter_fasta_records()`

- **Purpose:** Keep the FASTA records whose header matches a regular expression and write them to a new file.
- **Method:** Calls `genomepy.files.filter_fasta`, which returns the names of the removed sequences; the list is sorted. The output must differ from the input and must not exist.
- **API:** `dnakit.databases.filter_fasta_records(fasta_path, output_path, *, pattern=".*", invert_match=False, verbose=True)`.
- **Input:** Path to a FASTA file; a new output path; `pattern`; `invert_match`.
- **Output:** `output_path`, `pattern`, `invert_match`, `kept_count`, `removed` (sorted) and `provenance`.
- **Example code:**

```python
from pathlib import Path

from dnakit.databases import filter_fasta_records

Path("kept.fa").unlink(missing_ok=True)
result = filter_fasta_records("toy_n.fa", "kept.fa", pattern="chr", verbose=False)
print(result["removed"], result["kept_count"])
```

- **Example output:**

```text
['scaf_9'] 2
```

## 22) `dnakit.databases.sample_random_regions()`

- **Purpose:** Sample fixed-length regions from a FASTA file with a fixed random seed.
- **Method:** Calls `get_random_sequences` on a `genomepy.Genome`. The global Python random state is saved before, set from `seed`, and restored afterwards, so the caller's random numbers are unchanged. Regions whose N fraction exceeds `max_n` are rejected.
- **API:** `dnakit.databases.sample_random_regions(fasta_path, *, n=10, length=200, chroms=None, max_n=0.1, seed=0, verbose=True)`.
- **Input:** Path to a FASTA file; `n`, `length`, optional `chroms`, `max_n` (0 to 1) and `seed`.
- **Output:** `regions` (`chrom`, `start`, `end`; 0-based half-open), `region_count`, `seed` and `provenance`.
- **Example code:**

```python
import random
from pathlib import Path

from dnakit.databases import sample_random_regions

rng = random.Random(5)
chr1 = "".join(rng.choice("ACGT") for _ in range(3000))
Path("random.fa").write_text(f">chr1\n{chr1}\n>scaf_9\nNNNNACGTAC\n", encoding="utf-8")
result = sample_random_regions("random.fa", n=3, length=20, seed=7, verbose=False)
print(result["region_count"], result["seed"])
for region in result["regions"]:
    print(region["chrom"], region["start"], region["end"])
```

- **Example output:**

```text
3 7
chr1 215 235
chr1 1596 1616
chr1 1089 1109
```

## 23) `dnakit.databases.fetch_ncbi_assembly_metadata()`

- **Purpose:** Fetch NCBI Datasets assembly metadata for GCF or GCA accessions.
- **Method:** Calls `ncbi.datasets.metadata.genome.get_assembly_metadata_by_asm_accessions`, which needs network access. Warnings and errors that upstream prints are returned in `messages`; more than `max_records` raises `ConfigurationError`.
- **API:** `dnakit.databases.fetch_ncbi_assembly_metadata(accessions, *, max_records=100, verbose=True)`.
- **Input:** A list of 1 to 100 accessions, such as `GCF_000001405.40`; `max_records`.
- **Output:** `records`, `record_count`, `messages`, `requested` and `provenance`.
- **Example code (not run locally):**

```python
from dnakit.databases import fetch_ncbi_assembly_metadata

result = fetch_ncbi_assembly_metadata(["GCF_000001405.40"], max_records=5)
print(result["record_count"])
```

- **Example output (not run locally):**

```text
Not run locally: this function needs the NCBI Datasets network service; tests use a fake client.
```

## 24) `dnakit.databases.fetch_ncbi_assemblies_by_taxon()`

- **Purpose:** List NCBI Datasets assembly metadata for a taxonomy ID or species name.
- **Method:** Calls `ncbi.datasets.metadata.genome.get_assembly_metadata_by_taxon`, which needs network access. The output has the same shape as `fetch_ncbi_assembly_metadata`.
- **API:** `dnakit.databases.fetch_ncbi_assemblies_by_taxon(taxon, *, max_records=100, verbose=True)`.
- **Input:** An integer tax ID such as `9606`, or a name such as `Homo sapiens`; `max_records`.
- **Output:** `records`, `record_count`, `messages`, `requested` and `provenance`.
- **Example code (not run locally):**

```python
from dnakit.databases import fetch_ncbi_assemblies_by_taxon

result = fetch_ncbi_assemblies_by_taxon(9606, max_records=5)
print(result["record_count"])
```

- **Example output (not run locally):**

```text
Not run locally: this function needs the NCBI Datasets network service; tests use a fake client.
```

## 25) `dnakit.databases.fetch_ncbi_assemblies_by_bioproject()`

- **Purpose:** List NCBI Datasets assembly metadata for BioProject accessions.
- **Method:** Calls `ncbi.datasets.metadata.genome.get_assembly_metadata_by_bioproject_accessions`, which needs network access.
- **API:** `dnakit.databases.fetch_ncbi_assemblies_by_bioproject(bioproject_accessions, *, max_records=100, verbose=True)`.
- **Input:** A list of BioProject accessions starting with `PRJ` (such as `PRJNA168`); `max_records`.
- **Output:** `records`, `record_count`, `messages`, `requested` and `provenance`.
- **Example code (not run locally):**

```python
from dnakit.databases import fetch_ncbi_assemblies_by_bioproject

result = fetch_ncbi_assemblies_by_bioproject(["PRJNA168"], max_records=5)
print(result["record_count"])
```

- **Example output (not run locally):**

```text
Not run locally: this function needs the NCBI Datasets network service; tests use a fake client.
```

## 26) `dnakit.databases.select_ncbi_gene_fields()`

- **Purpose:** Select `gene_id` and `symbol` from one NCBI gene record.
- **Method:** Builds a `V1GeneDescriptor` from the mapping and calls `ncbi.datasets.metadata.gene.gene_values_by_fields`. Only these two fields exist, so other names raise an error. No network access.
- **API:** `dnakit.databases.select_ncbi_gene_fields(gene, fields, *, verbose=True)`.
- **Input:** A mapping with `gene_id` and `symbol`; `fields` is a non-empty subset of `gene_id` and `symbol`.
- **Output:** `fields`, `values` and `provenance`.
- **Example code:**

```python
from dnakit.databases import select_ncbi_gene_fields

gene = {"gene_id": "7157", "symbol": "TP53", "description": "tumor protein p53"}
print(select_ncbi_gene_fields(gene, ["gene_id", "symbol"], verbose=False)["values"])
```

- **Example output:**

```text
{'gene_id': '7157', 'symbol': 'TP53'}
```

## 27) `dnakit.databases.select_ncbi_assembly_fields()`

- **Purpose:** Select top-level fields from one assembly record.
- **Method:** Builds a `V1AssemblyDatasetDescriptor` and calls `ncbi.datasets.metadata.genome.assembly_values_by_fields`. An empty `fields` list returns every top-level field. No network access.
- **API:** `dnakit.databases.select_ncbi_assembly_fields(assembly, fields=(), *, verbose=True)`.
- **Input:** An assembly record mapping with string keys; `fields` is a list of top-level field names.
- **Output:** `fields`, `values`, `missing_fields` (requested names absent from the record) and `provenance`.
- **Example code:**

```python
from dnakit.databases import select_ncbi_assembly_fields

assembly = {"accession": "GCF_000001405.40", "organism": {"tax_id": 9606}}
print(select_ncbi_assembly_fields(assembly, ["accession", "organism"], verbose=False)["values"])
```

- **Example output:**

```text
{'accession': 'GCF_000001405.40', 'organism': {'tax_id': 9606}}
```

## How these functions are verified

- Test file: `others/tests/unit/test_full_G8.py`, 36 tests (including the parametrized check that every public function takes a keyword-only `verbose`).
- Each wrapper has at least one test that calls the real upstream package and compares the result with a hand calculation (FASTA lengths, interval merges, GTF column values, BED gap coordinates, Ensembl download paths) or with a direct upstream call.
- The NCBI Datasets assembly functions run the real upstream pagination generator with a fake API client; no test touches the network.
- Each module has a missing-package test that sets the package's `sys.modules` entry to `None` and checks for `BackendUnavailableError` with the extra name in the message.

## Current scope

- 27 functions are wrapped: 15 from pyensembl (FASTA reading, interval merging and nearest locus, identifier and alias normalization, version matching, reference names and Ensembl download URLs), 7 from genomepy (annotation read, write and filter, gap and sizes files, FASTA filtering, random region sampling), and 5 from NCBI Datasets (3 assembly metadata queries and 2 field selectors).
- Not wrapped: pyensembl's `Genome` and `EnsemblRelease` classes, cache and dated-release helpers and the CLI functions (they need downloads or a local cache directory); genomepy's `Genome` method family (`track2fasta` is wrapped as `write_bed_sequences`), the installer and provider queries (`install_genome` is wrapped as `install_genome_assembly`), functions that need UCSC tools, and configuration and plugin functions; NCBI Datasets' openapi model utilities, gRPC registration and `print_*` output functions.
- The examples import from `dnakit.databases`.
