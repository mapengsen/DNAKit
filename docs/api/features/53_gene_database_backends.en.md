# Optional gene-level database backends

This group of functions in `dnakit.databases` is a set of gene-level adapters over optional packages: [gget](https://github.com/pachterlab/gget), [pyensembl](https://github.com/openvax/pyensembl), [genomepy](https://github.com/vanheeringen-lab/genomepy), [biothings_client](https://github.com/biothings/biothings_client.py), [mygene](https://github.com/biothings/mygene.py), [myvariant](https://github.com/biothings/myvariant.py), [pybiomart](https://github.com/jrderuiter/pybiomart), and [ncbi-datasets-pylib](https://github.com/ncbi/datasets). Local functions read only the GTF, FASTA, BED, or VCF files passed explicitly, and compute gene loci, sequence slices, mutation windows, and HGVS identifiers; network functions call the Ensembl, MyGene.info, MyVariant.info, BioMart, or NCBI services and do no local computation. All dependencies are optional: a package is imported only when its function is called, and a missing package raises an error with a `pip install "dnakit[...]"` hint.

The local examples use only the four small files below. They are hand-written demonstration data, not a real genome.

## Example data

**Example file:** `genes.gtf`

```text
##description: tiny Ensembl-style annotation for DNAKit examples
chr1	test	gene	11	200	.	+	.	gene_id "ENSG00000000001"; gene_name "GENEA"; gene_biotype "protein_coding";
chr1	test	transcript	11	200	.	+	.	gene_id "ENSG00000000001"; gene_name "GENEA"; transcript_id "ENST00000000011"; transcript_name "GENEA-201"; transcript_biotype "protein_coding";
chr1	test	exon	11	60	.	+	.	gene_id "ENSG00000000001"; gene_name "GENEA"; transcript_id "ENST00000000011"; transcript_name "GENEA-201"; exon_number "1"; exon_id "ENSE00000000001";
chr1	test	exon	121	200	.	+	.	gene_id "ENSG00000000001"; gene_name "GENEA"; transcript_id "ENST00000000011"; transcript_name "GENEA-201"; exon_number "2"; exon_id "ENSE00000000002";
chr1	test	gene	301	380	.	-	.	gene_id "ENSG00000000002"; gene_name "GENEB"; gene_biotype "lncRNA";
chr1	test	transcript	301	380	.	-	.	gene_id "ENSG00000000002"; gene_name "GENEB"; transcript_id "ENST00000000022"; transcript_name "GENEB-201"; transcript_biotype "lncRNA";
chr1	test	exon	301	380	.	-	.	gene_id "ENSG00000000002"; gene_name "GENEB"; transcript_id "ENST00000000022"; transcript_name "GENEB-201"; exon_number "1"; exon_id "ENSE00000000003";
```

**Example file:** `genome.fa`

```text
>chr1 example contig
ACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGT
TGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCA
ACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGT
TGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCA
ACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGT
TGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCA
ACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCA
```

**Example file:** `regions.bed`

```text
chr1	10	20	r1	0	+
chr1	300	310	r2	0	-
```

**Example file:** `sample.vcf`

```text
##fileformat=VCFv4.2
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
1	10	rs1	A	G,T	.	PASS	.
chr1	20	.	C	.	.	.	.
chrUn	5	.	C	T	.	.	.
chr2	30	.	GA	G	.	.	.
```

The gene model has two artificial genes: GENEA is on the plus strand of `chr1:11-200`, with exons 11-60 and 121-200; GENEB is on the minus strand of `chr1:301-380`. The genome is a 400 bp periodic sequence, so slices can be checked by hand.

## 1) `dnakit.databases.search_ensembl_genes()`

- **Purpose:** Search Ensembl genes by free-text words.
- **Method:** Call gget's Ensembl search, request JSON output, and limit the number of returned rows.
- **API:** `dnakit.databases.search_ensembl_genes(searchwords, *, species="homo_sapiens", release=None, andor="or", max_results=100)`.
- **Input:** One word or a list of words; `species` in `genus_species` form; `andor` is `"or"` or `"and"`.
- **Output:** A dictionary with `count`, `results`, and `provenance`. Each result has the Ensembl ID, gene name, biotype, and description.
- **Example code (not run locally):**

```python
from dnakit.databases import search_ensembl_genes

result = search_ensembl_genes(["CDK2"], max_results=5)
print(result["count"], [row["gene_name"] for row in result["results"]])
```

- **Example output (not run locally):** This needs the live Ensembl service, and results change with Ensembl releases, so no fixed output is shown. `results` is a list whose items contain `ensembl_id` and `gene_name`.

## 2) `dnakit.databases.fetch_ensembl_gene_info()`

- **Purpose:** Fetch Ensembl gene metadata, including coordinates, strand, canonical transcript, and transcript lists.
- **Method:** Call gget's `info`. By default NCBI and UniProt are not queried, so only Ensembl data is returned and fewer requests are made.
- **API:** `dnakit.databases.fetch_ensembl_gene_info(ensembl_ids, *, external_sources=False)`.
- **Input:** A list of Ensembl IDs, at most 200. With `external_sources=True`, NCBI and UniProt are also queried.
- **Output:** A `records` dictionary keyed by the requested IDs. Each value is gget's metadata or `None` when the ID was not found. Also `requested` and `found` counts.
- **Example code (not run locally):**

```python
from dnakit.databases import fetch_ensembl_gene_info

result = fetch_ensembl_gene_info(["ENSG00000123374"])
print(result["records"]["ENSG00000123374"]["ensembl_gene_name"])
```

- **Example output (not run locally):** This needs the live Ensembl service. On success the gene name is `CDK2`; no fixed output is shown.

## 3) `dnakit.databases.fetch_ensembl_sequences()`

- **Purpose:** Fetch nucleotide sequences of genes or transcripts, or their protein translations, parsed into records.
- **Method:** Call gget's `seq`. gget returns alternating FASTA header and sequence lines; this function pairs them and checks the layout.
- **API:** `dnakit.databases.fetch_ensembl_sequences(ensembl_ids, *, translate=False, isoforms=False, max_records=100)`.
- **Input:** A list of Ensembl IDs. `translate=True` returns protein sequences; `isoforms=True` returns all transcripts.
- **Output:** A dictionary with `count` and `records`. Each record has `id`, `description`, and `sequence`.
- **Example code (not run locally):**

```python
from dnakit.databases import fetch_ensembl_sequences

result = fetch_ensembl_sequences(["ENSG00000123374"], translate=True)
print(result["records"][0]["id"], result["records"][0]["sequence"][:10])
```

- **Example output (not run locally):** This needs the live Ensembl or UniProt service, so no fixed output is shown. On success each record's `id` is an Ensembl ID and `sequence` is a protein sequence.

## 4) `dnakit.databases.apply_cdna_mutations()`

- **Purpose:** Apply one `c.`-notation mutation to each sequence and return the sequence window around the mutation.
- **Method:** Call gget's `mutate` once per sequence and mutation, so each output maps to its input. When the reference base at the position does not match, gget drops the mutation; this function reports it as `dropped`.
- **API:** `dnakit.databases.apply_cdna_mutations(sequences, mutations, *, flank=30)`.
- **Input:** `sequences` maps sequence IDs to DNA strings. `mutations` maps sequence IDs to one `c.` mutation, such as `"c.2C>T"`, with 1-based positions. `flank` is the number of bases kept on each side of the mutation.
- **Output:** A `mutations` list. Each item has `seq_id`, `mutation`, `status` (`applied` or `dropped`), and `mutant_window`. Also an `applied` count.
- **Example code:**

```python
from dnakit.databases import apply_cdna_mutations

result = apply_cdna_mutations({"s1": "ACGTTGCA" * 5}, {"s1": "c.2C>T"}, flank=3)
print(result["mutations"][0]["status"], result["mutations"][0]["mutant_window"])
```

- **Example output:**

```text
applied ATGTT
```

## 5) `dnakit.databases.index_gene_models()`

- **Purpose:** Index a local Ensembl-style GTF with pyensembl and count its contigs, genes, transcripts, and exons.
- **Method:** pyensembl parses the GTF into an SQLite index in the cache directory. The index name includes the first 16 hex digits of the GTF's SHA-256 digest, so an edited GTF gets a new index and never reuses a stale one.
- **API:** `dnakit.databases.index_gene_models(gtf, *, cache_dir=None)`.
- **Input:** The GTF path. `cache_dir` is the index cache directory; when omitted, pyensembl's default cache is used.
- **Output:** A dictionary with `contigs`, `genes`, `transcripts`, and `exons` counts, and `provenance`.
- **Example code:**

```python
from dnakit.databases import index_gene_models

summary = index_gene_models("genes.gtf", cache_dir="ensembl_cache")
print(summary["contigs"], summary["genes"], summary["transcripts"], summary["exons"])
```

- **Example output:**

```text
1 2 2 3
```

## 6) `dnakit.databases.query_genes_by_name()`

- **Purpose:** Find genes in a local GTF by exact gene name.
- **Method:** Call pyensembl's `gene_ids_of_gene_name`. pyensembl raises `ValueError` for a name with no match; this function converts that to an empty list.
- **API:** `dnakit.databases.query_genes_by_name(gtf, gene_names, *, cache_dir=None, max_results=100)`.
- **Input:** The GTF path, a list of gene names, and `max_results`, the upper limit on total matches.
- **Output:** A `results` dictionary keyed by the queried names. Each value is a list of gene records (`gene_id`, `contig`, `start`, `end`, `strand`, `biotype`, `transcript_ids`). Also a `matches` count.
- **Example code:**

```python
from dnakit.databases import query_genes_by_name

result = query_genes_by_name("genes.gtf", ["GENEA", "NOPE"], cache_dir="ensembl_cache")
print(result["results"]["GENEA"][0]["gene_id"], result["results"]["NOPE"])
```

- **Example output:**

```text
ENSG00000000001 []
```

## 7) `dnakit.databases.query_genes_at_locus()`

- **Purpose:** List genes that overlap a genomic locus.
- **Method:** Call pyensembl's `genes_at_locus`. Coordinates are 1-based and inclusive. When `end` is given, genes overlapping the whole interval are returned. `strand` can restrict the result to `"+"` or `"-"`.
- **API:** `dnakit.databases.query_genes_at_locus(gtf, contig, start, end=None, *, strand=None, cache_dir=None, max_results=100)`.
- **Input:** The GTF path, contig name, start, optional end, and optional strand. A contig that is not in the annotation raises an error.
- **Output:** A dictionary with `count` and `genes`. Each gene has the same record format as in `query_genes_by_name()`.
- **Example code:**

```python
from dnakit.databases import query_genes_at_locus

result = query_genes_at_locus("genes.gtf", "chr1", 190, 310, cache_dir="ensembl_cache")
print([gene["gene_name"] for gene in result["genes"]])
```

- **Example output:**

```text
['GENEA', 'GENEB']
```

## 8) `dnakit.databases.fetch_gene_sequence()`

- **Purpose:** Extract a gene's genomic locus from a local genome FASTA, oriented by the gene's strand.
- **Method:** pyensembl finds the gene coordinates, then `Genome.sequence()` reads the sequence. Minus-strand genes return the reverse complement. A sequence longer than `max_length` raises an error.
- **API:** `dnakit.databases.fetch_gene_sequence(gtf, genome_fasta, gene_id, *, cache_dir=None, max_length=1000000)`.
- **Input:** The GTF path, the genome FASTA path, an Ensembl gene ID, and the length limit.
- **Output:** A dictionary with `gene_id`, `gene_name`, `contig`, `start`, `end`, `strand`, `length`, and `sequence`.
- **Example code:**

```python
from dnakit.databases import fetch_gene_sequence

result = fetch_gene_sequence("genes.gtf", "genome.fa", "ENSG00000000002", cache_dir="ensembl_cache")
print(result["length"], result["sequence"][:20])
```

- **Example output:**

```text
80 ACGTTGCAACGTTGCAACGT
```

## 9) `dnakit.databases.read_genome_summary()`

- **Purpose:** Summarize the contigs, lengths, and N-base counts of a local genome FASTA.
- **Method:** Call genomepy's `Genome`. It reads a pyfaidx index and writes `.fai`, `.fa.sizes`, and `.gaps.bed` files next to the FASTA. The FASTA itself is not modified. `na` metadata from the README is returned as `None`.
- **API:** `dnakit.databases.read_genome_summary(fasta, *, max_contigs=100000)`.
- **Input:** The FASTA path. `max_contigs` is the record-count limit.
- **Output:** A dictionary with `name`, `tax_id`, `assembly_accession`, `contig_count`, `total_length`, and `contigs`.
- **Example code:**

```python
from dnakit.databases import read_genome_summary

summary = read_genome_summary("genome.fa")
print(summary["name"], summary["contig_count"], summary["total_length"])
```

- **Example output:**

```text
genome 1 400
```

## 10) `dnakit.databases.fetch_spliced_sequence()`

- **Purpose:** Join intervals from one contig in order, as for spliced exons.
- **Method:** Call genomepy's `get_spliced_seq`, which is inherited from pyfaidx. Intervals are 1-based and inclusive, so `[11, 60]` covers bases 11 to 60. With `reverse_complement=True`, the joined sequence is reverse-complemented.
- **API:** `dnakit.databases.fetch_spliced_sequence(fasta, contig, intervals, *, reverse_complement=False, max_length=1000000)`.
- **Input:** The FASTA path, contig name, a list of `[start, end]` intervals, and the length limit. An interval past the end of the contig raises an error.
- **Output:** A dictionary with `contig`, `intervals`, `length`, and `sequence`.
- **Example code:**

```python
from dnakit.databases import fetch_spliced_sequence

result = fetch_spliced_sequence("genome.fa", "chr1", [[11, 60], [121, 200]])
print(result["length"], result["sequence"][:20])
```

- **Example output:**

```text
130 GTTGCAACGTTGCAACGTTG
```

## 11) `dnakit.databases.write_bed_sequences()`

- **Purpose:** Extract sequences for BED regions from a genome FASTA and write them to a new FASTA file.
- **Method:** Call genomepy's `track2fasta`. BED coordinates are 0-based and half-open. With `stranded=True`, minus-strand regions are returned as reverse complements. If `output` already exists, the function raises an error and does not overwrite it.
- **API:** `dnakit.databases.write_bed_sequences(fasta, bed, output, *, stranded=False, extend_up=0, extend_down=0, max_regions=100000)`.
- **Input:** The FASTA path, the BED path, and the output path. `extend_up` and `extend_down` extend each region by that many bases in the strand direction.
- **Output:** A dictionary with `regions`, `records`, `skipped`, `output`, and `provenance`.
- **Example code:**

```python
from pathlib import Path
from dnakit.databases import write_bed_sequences

output = Path("regions.fa")
output.unlink(missing_ok=True)  # the example is re-runnable; the function itself never overwrites
result = write_bed_sequences("genome.fa", "regions.bed", output, stranded=True)
print(result["records"], output.read_text().splitlines()[1])
```

- **Example output:**

```text
2 GTTGCAACGT
```

## 12) `dnakit.databases.install_genome_assembly()`

- **Purpose:** Download and install a genome assembly through genomepy, optionally with its annotation.
- **Method:** Call `genomepy.install_genome`. It downloads large files from the chosen provider, so it takes time and disk space.
- **API:** `dnakit.databases.install_genome_assembly(name, *, provider=None, genomes_dir=None, annotation=False, threads=1)`.
- **Input:** An assembly name such as `hg38`, a provider such as `ucsc` or `ensembl`, the installation directory, whether to download the annotation, and the number of threads (1 to 64).
- **Output:** A dictionary with `name`, `provider`, `fasta`, and `annotation_gtf`.
- **Example code (not run locally):**

```python
from dnakit.databases import install_genome_assembly

result = install_genome_assembly("hg38", provider="ucsc", annotation=True)
print(result["fasta"])
```

- **Example output (not run locally):** This downloads several gigabytes of genome data, so it is not run here. On success `fasta` points to the FASTA file in the installation directory.

## 13) `dnakit.databases.format_hgvs_id()`

- **Purpose:** Format one VCF-style variant as a genomic HGVS identifier.
- **Method:** Call biothings_client's variant client `format_hgvs`. This is pure string formatting and makes no network request. The identifier does not carry a genome build, so use the matching `assembly` in later MyVariant queries.
- **API:** `dnakit.databases.format_hgvs_id(chrom, pos, ref, alt)`.
- **Input:** A human chromosome name, with or without the `chr` prefix; a 1-based position; the reference and alternate bases, limited to A, C, G, T, and N.
- **Output:** A dictionary with `hgvs_id`, `chrom`, `pos`, `ref`, and `alt`.
- **Example code:**

```python
from dnakit.databases import format_hgvs_id

print(format_hgvs_id("chr9", 107620835, "G", "A")["hgvs_id"])
```

- **Example output:**

```text
chr9:g.107620835G>A
```

## 14) `dnakit.databases.hgvs_ids_from_vcf()`

- **Purpose:** Produce a genomic HGVS identifier for each sequence alternate allele in a VCF.
- **Method:** Parse the VCF line by line and pass each sequence alternate allele to `format_hgvs`. Records are counted in `skipped` when the contig is not a human chromosome, when an alternate allele is `.` or structural, or when it equals the reference. The upstream VCF iterator does not perform these checks, so it is not used directly.
- **API:** `dnakit.databases.hgvs_ids_from_vcf(path, *, max_records=100000)`.
- **Input:** The VCF path and the record limit.
- **Output:** A dictionary with `records`, `hgvs_ids`, and `skipped`.
- **Example code:**

```python
from dnakit.databases import hgvs_ids_from_vcf

result = hgvs_ids_from_vcf("sample.vcf")
print(result["hgvs_ids"])
print(result["skipped"])
```

- **Example output:**

```text
['chr1:g.10A>G', 'chr1:g.10A>T', 'chr2:g.31del']
{'unsupported_contig': 1, 'non_sequence_allele': 1, 'same_as_ref': 0}
```

## 15) `dnakit.databases.query_mygene_genes()`

- **Purpose:** Search MyGene.info gene records with a query string.
- **Method:** Call mygene's `query`, and return the total count and the hits, limited to `size` rows.
- **API:** `dnakit.databases.query_mygene_genes(query, *, species="human", fields="symbol,name,entrezgene,ensembl.gene,taxid", size=10)`.
- **Input:** A MyGene query string, such as `"symbol:CDK2"`, plus the species, field list, and number of rows (1 to 1000).
- **Output:** A dictionary with `total` and `hits`.
- **Example code (not run locally):**

```python
from dnakit.databases import query_mygene_genes

result = query_mygene_genes("symbol:CDK2", size=1)
print(result["hits"][0]["symbol"])
```

- **Example output (not run locally):** This needs the live MyGene.info service. On success the first hit's `symbol` is `CDK2`.

## 16) `dnakit.databases.fetch_mygene_gene()`

- **Purpose:** Fetch one MyGene.info gene record by Entrez or Ensembl ID.
- **Method:** Call mygene's `getgene`. A missing gene returns `found=False` instead of raising an error.
- **API:** `dnakit.databases.fetch_mygene_gene(gene_id, *, species="human", fields="symbol,name,entrezgene,ensembl.gene,taxid")`.
- **Input:** A gene ID such as `"1017"`, and the species.
- **Output:** A dictionary with `gene_id`, `found`, and `gene`.
- **Example code (not run locally):**

```python
from dnakit.databases import fetch_mygene_gene

result = fetch_mygene_gene("1017")
print(result["found"], result["gene"]["symbol"])
```

- **Example output (not run locally):** This needs the live MyGene.info service. On success it prints `True` and `CDK2`.

## 17) `dnakit.databases.map_gene_identifiers()`

- **Purpose:** Map a batch of identifiers of one type, such as Entrez IDs, to MyGene.info gene records.
- **Method:** Call mygene's `querymany`. Each result keeps a `query` field so it can be matched to its input.
- **API:** `dnakit.databases.map_gene_identifiers(identifiers, *, scopes="entrezgene", species="human", fields="symbol,name,taxid")`.
- **Input:** A list of identifiers (at most 1000), the identifier type in `scopes`, and the species.
- **Output:** A dictionary with `requested` and `results`.
- **Example code (not run locally):**

```python
from dnakit.databases import map_gene_identifiers

result = map_gene_identifiers(["1017", "7157"])
print([row["symbol"] for row in result["results"]])
```

- **Example output (not run locally):** This needs the live MyGene.info service. On success the symbols include `CDK2` and `TP53`.

## 18) `dnakit.databases.fetch_variant_annotations()`

- **Purpose:** Fetch MyVariant.info annotations for a list of genomic HGVS identifiers.
- **Method:** Call myvariant's `getvariants`, and pass `assembly` to the service as a request parameter.
- **API:** `dnakit.databases.fetch_variant_annotations(hgvs_ids, *, assembly="hg38", fields="dbsnp.rsid,clinvar.rcv.clinical_significance,cadd.phred,dbnsfp.genename")`.
- **Input:** A list of HGVS identifiers (at most 1000) such as `"chr9:g.107620835G>A"`. `assembly` is `"hg19"` or `"hg38"`.
- **Output:** A dictionary with `requested`, `results`, and `provenance`. Identifiers that are not found carry a `notfound` flag in `results`.
- **Example code (not run locally):**

```python
from dnakit.databases import fetch_variant_annotations

result = fetch_variant_annotations(["chr9:g.107620835G>A"], assembly="hg19")
print(result["results"][0]["_id"])
```

- **Example output (not run locally):** This needs the live MyVariant.info service. On success it prints `chr9:g.107620835G>A`.

## 19) `dnakit.databases.query_variant_annotations()`

- **Purpose:** Search MyVariant.info variant records with a query string.
- **Method:** Call myvariant's `query`, and return the total count and the hits, limited to `size` rows.
- **API:** `dnakit.databases.query_variant_annotations(query, *, assembly="hg38", fields=..., size=10)`.
- **Input:** A query string such as `"dbsnp.rsid:rs58991260"`, the `assembly`, and the number of rows (1 to 1000).
- **Output:** A dictionary with `total` and `hits`.
- **Example code (not run locally):**

```python
from dnakit.databases import query_variant_annotations

result = query_variant_annotations("dbsnp.rsid:rs58991260", size=1)
print(result["total"])
```

- **Example output (not run locally):** This needs the live MyVariant.info service. The number of hits changes as the database is updated, so no fixed output is shown.

## 20) `dnakit.databases.query_biomart_dataset()`

- **Purpose:** Query a BioMart dataset for the requested attributes and filters.
- **Method:** Call pybiomart's `Dataset.query`, which builds the BioMart XML request and parses the TSV response. Empty cells become `None`. Importing pybiomart creates a `.pybiomart.sqlite` request cache in the current working directory.
- **API:** `dnakit.databases.query_biomart_dataset(dataset, attributes, *, filters=None, host="http://www.ensembl.org", only_unique=True, use_attr_names=False, max_rows=10000)`.
- **Input:** A dataset name such as `hsapiens_gene_ensembl`, a list of attributes, an optional filter mapping, the service host, and the row limit. The service returns the full result before `max_rows` is checked.
- **Output:** A dictionary with `columns`, `row_count`, and `rows`.
- **Example code (not run locally):**

```python
from dnakit.databases import query_biomart_dataset

result = query_biomart_dataset(
    "hsapiens_gene_ensembl",
    ["ensembl_gene_id", "external_gene_name"],
    filters={"chromosome_name": "21"},
    max_rows=5,
)
print(result["rows"][0])
```

- **Example output (not run locally):** This needs the live Ensembl BioMart service. On success each item of `rows` is a dictionary with the two requested columns.

## 21) `dnakit.databases.list_biomart_attributes()`

- **Purpose:** List the attributes that a BioMart dataset exposes.
- **Method:** Call pybiomart's `Dataset.list_attributes`, which returns the attribute name, display name, and description.
- **API:** `dnakit.databases.list_biomart_attributes(dataset, *, host="http://www.ensembl.org", max_rows=5000)`.
- **Input:** A dataset name, the service host, and the row limit.
- **Output:** A dictionary with `columns`, `count`, and `attributes`.
- **Example code (not run locally):**

```python
from dnakit.databases import list_biomart_attributes

result = list_biomart_attributes("hsapiens_gene_ensembl")
print(result["count"])
```

- **Example output (not run locally):** This needs the live Ensembl BioMart service. The attribute count changes between Ensembl releases, so no fixed output is shown.

## 22) `dnakit.databases.list_biomart_datasets()`

- **Purpose:** List the datasets in one BioMart mart.
- **Method:** Use pybiomart's `Server.marts` to find the mart, then call `Mart.list_datasets`. A mart name that the service does not list raises an error.
- **API:** `dnakit.databases.list_biomart_datasets(mart="ENSEMBL_MART_ENSEMBL", *, host="http://www.ensembl.org", max_rows=1000)`.
- **Input:** The mart name, the service host, and the row limit.
- **Output:** A dictionary with `mart`, `columns`, `count`, and `datasets`.
- **Example code (not run locally):**

```python
from dnakit.databases import list_biomart_datasets

result = list_biomart_datasets(max_rows=5)
print(result["datasets"][0]["name"])
```

- **Example output (not run locally):** This needs the live Ensembl BioMart service. The dataset list changes between releases, so no fixed output is shown.

## 23) `dnakit.databases.fetch_ncbi_gene_by_id()`

- **Purpose:** Fetch NCBI gene metadata for NCBI Gene IDs.
- **Method:** Call ncbi-datasets-pylib's `GeneApi.gene_metadata_by_id`, with the server-side `limit` set to `max_genes`. The openapi response model is converted to plain data with `to_dict()`.
- **API:** `dnakit.databases.fetch_ncbi_gene_by_id(gene_ids, *, max_genes=100)`.
- **Input:** A list of positive integer Gene IDs, at most 200; `max_genes` is the result limit (1 to 1000).
- **Output:** A dictionary with `requested`, `total_count`, and `genes`. Each item of `genes` is an NCBI gene match record.
- **Example code (not run locally):**

```python
from dnakit.databases import fetch_ncbi_gene_by_id

result = fetch_ncbi_gene_by_id([1017], max_genes=1)
print(result["genes"][0]["gene"]["symbol"])
```

- **Example output (not run locally):** This needs the live NCBI Datasets service and the `ncbi-datasets-pylib` package. On success it prints the symbol `CDK2`.

## 24) `dnakit.databases.fetch_ncbi_gene_by_symbol()`

- **Purpose:** Fetch NCBI gene metadata by gene symbol and taxon.
- **Method:** Call `GeneApi.gene_metadata_by_tax_and_symbol`. Symbols and the taxon are format-checked first.
- **API:** `dnakit.databases.fetch_ncbi_gene_by_symbol(symbols, taxon, *, max_genes=100)`.
- **Input:** A list of gene symbols, a taxon such as `"human"` or `"9606"`, and the result limit.
- **Output:** A dictionary with `symbols`, `taxon`, `total_count`, and `genes`.
- **Example code (not run locally):**

```python
from dnakit.databases import fetch_ncbi_gene_by_symbol

result = fetch_ncbi_gene_by_symbol(["CDK2"], "human", max_genes=1)
print(result["genes"][0]["gene"]["gene_id"])
```

- **Example output (not run locally):** This needs the live NCBI Datasets service and the `ncbi-datasets-pylib` package. On success it prints the Gene ID for `CDK2`.

## 25) `dnakit.databases.fetch_ncbi_gene_orthologs()`

- **Purpose:** Fetch the NCBI ortholog set for one Gene ID.
- **Method:** Call `GeneApi.gene_orthologs_by_id`. When `taxon_filter` is omitted, no taxon filter is applied.
- **API:** `dnakit.databases.fetch_ncbi_gene_orthologs(gene_id, *, taxon_filter=None)`.
- **Input:** A positive integer Gene ID and an optional list of taxon IDs, such as `["9606"]`.
- **Output:** A dictionary with `gene_id`, `taxon_filter`, `ortholog_set_id`, `method`, and `genes`.
- **Example code (not run locally):**

```python
from dnakit.databases import fetch_ncbi_gene_orthologs

result = fetch_ncbi_gene_orthologs(1017, taxon_filter=["9606"])
print(result["method"], len(result["genes"]))
```

- **Example output (not run locally):** This needs the live NCBI Datasets service and the `ncbi-datasets-pylib` package. The output is the orthology method name and the number of orthologs, so no fixed output is shown.

## How these functions are verified

- The test file is `others/tests/unit/test_database_gene_backends.py`, with 37 tests: 33 pass and 4 are skipped. The skipped tests are two pybiomart tests (dependency conflict, see below), one ncbi-datasets-pylib signature test (package not installed), and one live test that needs `DNAKIT_LIVE_TESTS=1`.
- Local tests run the real upstream packages. pyensembl indexes a GTF in a temporary directory and its counts are compared with a direct `pyensembl.Genome` call. genomepy's spliced sequence is compared with direct `get_seq` calls. gget.mutate output is compared with a direct call. biothings_client HGVS strings are compared with the upstream `format_hgvs` and `get_hgvs_from_vcf`.
- Sequences, exon joins, mutation windows, and HGVS strings are also compared with hand-calculated slices. For example, GENEA is 190 bp, and GENEB equals the reverse complement of its locus.
- Network functions replace only the upstream client call and test parsing, missing values, and argument validation without network access. The ncbi functions are tested with a fake `ncbi.datasets` module that checks arguments, the server-side `limit`, and error conversion. The real-signature test runs when the package can be imported.
- A missing optional package raises `BackendUnavailableError` with the matching `bio-*` extra name.
- Setting `DNAKIT_LIVE_TESTS=1` runs one live test that checks CDK2 through Ensembl and MyGene.
- Input validation covers booleans, out-of-range numbers, malformed identifiers, overlong sequences, existing output files, and unknown contigs.

## Current scope

- Covered: gget's Ensembl search, gene information, sequences, and cDNA mutation windows; pyensembl's local GTF index, gene name and locus queries, and gene sequences; genomepy's FASTA summary, spliced sequences, BED sequence output, and assembly installation; biothings_client's HGVS formatting; MyGene and MyVariant queries; BioMart attributes, datasets, and queries; and NCBI Datasets gene metadata, symbol lookup, and orthologs.
- Not wrapped: gget's `ref` (it only returns FTP links), `blat`, `blast`, `cbio`, `cellxgene`, and other online analysis interfaces; pyensembl's `EnsemblRelease` download workflow, since this module uses only explicitly given local GTF and FASTA files; and biothings_client's gene and variant clients, which are the same services as MyGene and MyVariant, so only its HGVS tools are wrapped.
- Extras: `bio-gget` = `gget>=0.30,<1`; `bio-pyensembl` = `pyensembl>=2.27,<3`; `bio-genomepy` = `genomepy>=0.16,<1`; `bio-biothings` = `biothings_client>=0.5,<1`; `bio-mygene` = `mygene>=3.2,<4`; `bio-myvariant` = `myvariant>=1.0,<2`; `bio-pybiomart` = `pybiomart>=0.2,<1`; `bio-ncbi-datasets` = `ncbi-datasets-pylib>=16.6,<17`.
- ncbi-datasets-pylib (16.6.1) cannot be installed in this venv: PyPI publishes only an sdist, and its build fails because `pkg_resources` is missing. The three NCBI functions were written from that sdist's source (the `GeneApi` method signatures and the `limit` parameter) and have not been run against the real package.
- pybiomart cannot be imported in this venv: `cattrs` 26.2.1 needs a newer `attrs`, and the venv has attrs 21.4.0. The tests for the three BioMart functions are skipped in this environment.
- Live-service results change over time, so the examples that call live services are not run locally and no output is shown for them.
