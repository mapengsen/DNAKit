# H4 family-variant check: polars-bio, pyensembl, gtars, bioframe and five more packages

This page records the review of the `family-variant` ledger rows in group H4. Each row had been marked as a variant of an already wrapped function (an alias, a re-export, or a branch chosen by a parameter). Every claim was checked: aliases and re-exports became `duplicate-of`; variants that an existing wrapper really reaches through a parameter stay `family-variant`; variants no existing wrapper reaches got a new wrapper. The new wrappers are in `src/dnakit/annotation/gtars_family_functions.py`, `src/dnakit/annotation/bioframe_family_functions.py` and `src/dnakit/databases/pyensembl_family_functions.py`. They are exported from the area packages, and the examples import them from there.

- 132 rows were checked, across polars_bio, bioframe, gtars, pybedtools, gffutils, gtfparse and genomicranges (annotation area) and pyensembl and genomepy (databases area).
- 72 rows became `duplicate-of`, 48 stay `family-variant` after the parameter was verified, 8 got new wrappers, 2 are `internal` and 2 are `broken-upstream`.
- The tests are in `others/tests/unit/test_family_H4.py`. A kept `family-variant` row is checked by calling the existing wrapper with the claimed parameter and comparing with the upstream call. A `duplicate-of` row is checked by object identity or by reading the upstream source.
- The example outputs come from running the code, so the documentation checker should report 0 mismatches.

## Review results

| Package | Family-variant before | duplicate-of | family-variant (verified) | wrapped (new) | Other status |
| --- | --- | --- | --- | --- | --- |
| polars_bio | 52 | 22 | 28 | 0 | internal 2 |
| pyensembl | 42 | 38 | 1 | 3 | none |
| gffutils | 10 | 0 | 10 | 0 | none |
| genomepy | 9 | 5 | 4 | 0 | none |
| bioframe | 6 | 3 | 0 | 3 | none |
| gtars | 6 | 1 | 3 | 2 | none |
| pybedtools | 4 | 0 | 2 | 0 | broken-upstream 2 |
| gtfparse | 2 | 2 | 0 | 0 | none |
| genomicranges | 1 | 1 | 0 | 0 | none |

Main findings:

- 38 pyensembl rows are re-exports: paths such as `pyensembl.genome.normalize_chromosome` are the same object as the function in its defining module (checked with `is`), not variants. Only `find_species_by_reference`, defined in `pyensembl.reference_name`, stays as `ensembl_species_record(by='reference')`.
- `make_gtf_filename` and `make_fasta_filename` return the file name alone. No parameter of the URL wrappers returns only the file name, so they were wrapped.
- The claim for `is_canonical_source` was wrong. It is true for any file name in the official `dna/` directory, while the `canonical` field of `ensembl_genome_fasta_source_options` is true only for toplevel or primary_assembly DNA file names. For `.../dna/Homo_sapiens.GRCh38.cdna.all.fa.gz` the first is True and the second False.
- gtars `digest_sequence` and `load_fasta` return the alphabet gtars detects (`dna2bit`, `dna3bit`, `dnaio` and so on). The existing `sequence_digest` and `digest_fasta_sequences` do not, so both were wrapped.
- bioframe `select_mask`, `select_indices` and `select_labels` return row positions, not interval records, and `select_bioframe_region` does not carry them. One wrapper was added. The three aliases under `bioframe.ops` are the same objects as the top-level functions.
- pybedtools `random_intersection` and `random_intersection_bp` are not variants and do not work here: upstream `close_or_delete` throws `StopIteration` into an already closed generator, so every call fails. They are marked `broken-upstream`.
- polars_bio `scan_*` and `sink_*` are not reachable by a parameter: the source of `read_*` is `scan_*(...).collect()`, and `write_*` and `sink_*` end in the same writer. They became `duplicate-of`. `print_metadata_json` and `print_metadata_summary` only print the content of `get_metadata` and return None, so they are `internal`.
- gtfparse `parse_gtf_pandas` is not a pandas-returning variant: `parse_gtf` already returns pandas, and `parse_gtf_pandas` only adds progress callbacks. It became `duplicate-of`. genomicranges `parse_gtf` is the first step of `read_gtf` and became `duplicate-of` too.
- All 10 gffutils `merge_criteria` functions are reached by the `criteria` and `threshold` parameters of `merge_gffutils_features`. The test compares each with a direct `db.merge` call.
- genomepy `region_to_seq`, `regions_to_seqs`, `bed_to_seq` and `bed_to_seqs` are called by `Genome.track2fasta`. The `bed` parameter of `write_bed_sequences` reaches them with a regions file or a BED file.
- polars_bio `write_bam` and `write_cram` stay as `write_alignment_file(output_format=...)`. The upstream `sort_on_write` and `tag_type_overrides` arguments are not exposed.

## Shared conventions

- The `verbose` argument of every public function is keyword-only and defaults to `True`.
- Return values are JSON-serializable dicts with `status`, `operation`, the result fields and `provenance`.
- `pyensembl`, `gtars` and `bioframe` are imported only when a function runs. A missing dependency raises `BackendUnavailableError`.
- Invalid input raises `ConfigurationError`.
- The three pyensembl functions work on strings only and make no network request.

## 1) `ensembl_gtf_filename`

- **Purpose:** Return the file name of the gene annotation GTF for a species and Ensembl release.
- **Method:** Calls the upstream `ensembl_url_templates.make_gtf_filename`. The name equals the last segment of the URL from `ensembl_gtf_url`, but the URL wrapper has no parameter that returns the file name alone.
- **API:** `ensembl_gtf_filename(ensembl_release[required], species[required], *, verbose[optional])`
- **Input:** `ensembl_release` is an integer from 1 to 1000. `species` is a species name such as `homo_sapiens`.
- **Output:** A dict with `species`, `ensembl_release`, `filename` and `provenance`.

**Example code:**

```python
from dnakit.databases import ensembl_gtf_filename

print(ensembl_gtf_filename(75, "homo_sapiens", verbose=False)["filename"])
print(ensembl_gtf_filename(109, "homo_sapiens", verbose=False)["filename"])
```

**Example output:**

```text
Homo_sapiens.GRCh37.75.gtf.gz
Homo_sapiens.GRCh38.109.chr_patch_hapl_scaff.gtf.gz
```

## 2) `ensembl_sequence_fasta_filename`

- **Purpose:** Return the file name of the cDNA, ncRNA or protein FASTA for a species and release.
- **Method:** Calls the upstream `ensembl_url_templates.make_fasta_filename`. Upstream decides between the new directory layout (releases above 75) and the old one, which name files differently.
- **API:** `ensembl_sequence_fasta_filename(ensembl_release[required], species[required], sequence_type[required], *, verbose[optional])`
- **Input:** `sequence_type` is `cdna`, `ncrna` or `pep`. The other arguments are as above.
- **Output:** A dict with `species`, `ensembl_release`, `sequence_type`, `filename` and `provenance`.

**Example code:**

```python
from dnakit.databases import ensembl_sequence_fasta_filename

for kind in ("cdna", "ncrna", "pep"):
    print(kind, ensembl_sequence_fasta_filename(109, "homo_sapiens", kind, verbose=False)["filename"])
```

**Example output:**

```text
cdna Homo_sapiens.GRCh38.cdna.all.fa.gz
ncrna Homo_sapiens.GRCh38.ncrna.fa.gz
pep Homo_sapiens.GRCh38.pep.all.fa.gz
```

## 3) `check_ensembl_canonical_source`

- **Purpose:** Tell whether a URL is in the official Ensembl reference FASTA layout `.../fasta/<species>/dna/<file>`.
- **Method:** Calls the upstream `genome_fasta_cache.is_canonical_source`, which checks the scheme (https), the host and the path layout, not the file name. To read the DNA type and mask, use `ensembl_genome_fasta_source_options`.
- **API:** `check_ensembl_canonical_source(source[required], *, verbose[optional])`
- **Input:** `source` is a URL string.
- **Output:** A dict with `source`, `official_source` (boolean) and `provenance`.

**Example code:**

```python
from dnakit.databases import ensembl_genome_fasta_source_options
from dnakit.databases import check_ensembl_canonical_source

base = "https://ftp.ensembl.org/pub/release-109/fasta/homo_sapiens/dna/"
for name in ("Homo_sapiens.GRCh38.dna.toplevel.fa.gz", "Homo_sapiens.GRCh38.cdna.all.fa.gz"):
    url = base + name
    official = check_ensembl_canonical_source(url, verbose=False)["official_source"]
    dna = ensembl_genome_fasta_source_options(url, verbose=False)["canonical"]
    print(name, official, dna)
```

**Example output:**

```text
Homo_sapiens.GRCh38.dna.toplevel.fa.gz True True
Homo_sapiens.GRCh38.cdna.all.fa.gz True False
```

## 4) `sequence_digest_record`

- **Purpose:** Return the GA4GH sha512t24u and MD5 digests, the length and the gtars-detected alphabet of one sequence.
- **Method:** The sequence is uppercased and checked for DNA symbols, then the upstream `gtars.refget.digest_sequence` is called. `name` defaults to an empty string.
- **API:** `sequence_digest_record(sequence[required], *, name[optional], description[optional], verbose[optional])`
- **Input:** `sequence` is a non-empty DNA string (IUPAC symbols allowed). `name` and `description` are optional text.
- **Output:** A dict with `name`, `description`, `length`, `alphabet`, `sha512t24u`, `md5` and `provenance`.

**Example code:**

```python
from dnakit.annotation import sequence_digest_record

result = sequence_digest_record("acgtacgtnn", name="s1", verbose=False)
print(result["length"], result["alphabet"], result["sha512t24u"], result["md5"])
print(sequence_digest_record("ACGT", verbose=False)["alphabet"])
```

**Example output:**

```text
10 dna3bit 5rfV7JOOAEVoyCNHpPwNcJJwGyVAjmH5 a0e5e77e1f26e1b3d8a27baa3db84c40
dna2bit
```

## 5) `load_fasta_sequence_collection`

- **Purpose:** Load a FASTA file with gtars and return each sequence with its digests, alphabet and, optionally, its residues.
- **Method:** Calls the upstream `gtars.refget.load_fasta`. The collection digest equals the one from `digest_fasta_sequences`. Results are sorted by sequence name. A total above `max_bases` bases raises ConfigurationError.
- **API:** `load_fasta_sequence_collection(fasta[required], *, include_sequences[optional], max_rows[optional], max_bases[optional], verbose[optional])`
- **Input:** `fasta` is a local FASTA path. `include_sequences` defaults to `True`, `max_rows` to 100000 and `max_bases` to 100000000.
- **Output:** A dict with `sequence_count`, `total_bases`, `collection_digest`, `sequences` and `provenance`. Each sequence has `name`, `description`, `length`, `alphabet`, `sha512t24u`, `md5` and `sequence` (uppercase, no line breaks).

**Example file:** `t.fa`

```text
>a first record
ACGTAC
GGTT
>b
NNACGU
```

**Example code:**

```python
from dnakit.annotation import load_fasta_sequence_collection

result = load_fasta_sequence_collection("t.fa", verbose=False)
print(result["sequence_count"], result["total_bases"], result["collection_digest"])
for row in result["sequences"]:
    print(row["name"], row["length"], row["alphabet"], row["sequence"])
```

**Example output:**

```text
2 16 KvfiMp1LuIdD2NYCZudr72JOMy2jnJNr
a 10 dna2bit ACGTACGGTT
b 6 dnaio NNACGU
```

## 6) `select_bioframe_overlap_rows`

- **Purpose:** Return which rows of a BED file overlap a region, as a boolean mask, row positions and index labels.
- **Method:** Reads a BED3 to BED6 file and calls bioframe `select_mask`, `select_indices` and `select_labels` in turn. Intervals are 0-based half-open. A bare chromosome name selects every interval on it. The file is read with a default integer index, so `labels` equals `indices`.
- **API:** `select_bioframe_overlap_rows(source[required], region[required], *, max_rows[optional], verbose[optional])`
- **Input:** `source` is a local BED path. `region` is a UCSC-style region string. `max_rows` defaults to 100000.
- **Output:** A dict with `region`, `row_count`, `selected_count`, `mask` (one boolean per data row, in file order), `indices`, `labels` and `provenance`.

**Example file:** `a.bed`

```text
chr1	10	30
chr1	50	80
chr2	5	15
chr1	100	120
```

**Example code:**

```python
from dnakit.annotation import select_bioframe_overlap_rows

result = select_bioframe_overlap_rows("a.bed", "chr1:25-55", verbose=False)
print(result["mask"], result["indices"], result["labels"])
print(select_bioframe_overlap_rows("a.bed", "chr1", verbose=False)["indices"])
```

**Example output:**

```text
[True, True, False, False] [0, 1] [0, 1]
[0, 1, 3]
```
