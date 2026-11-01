# gtars transcript stores: coordinate conversion and lookup

This page covers the gtars `reftx` rows of work-list pack P24. A transcript store is one binary file built from transcript records, cdot JSON files and a MANE summary. The three coordinate rows (c. to g., n. to g., gene-based c. to g.) are one function with `coordinate_type` or `gene`; the lookup, chromosome-accession and strand rows are one lookup function; transcript adding, chromosome mapping, cdot loading and MANE loading are the parameters of one builder function. gtars 0.10.0 was tested.

## Example files

- **Example file:** `c.json`

```text
{"transcripts": {"NM_9.1": {"id": "NM_9.1", "gene_name": "GX", "contig": "NC_000001.11", "strand": 1, "cds_start": 150, "cds_end": 950, "exons": [[100, 400], [600, 1000]]}}}
```

- **Example file:** `mane.tsv`

```text
#NCBI_GeneID	Ensembl_Gene	HGNC_ID	symbol	name	RefSeq_nuc	RefSeq_prot	Ensembl_nuc	Ensembl_prot	MANE_status	GRCh38_chr	chr_start	chr_end	chr_strand
1	ENSG1	HGNC:1	GX	gene x	NM_9.1	NP_9.1	ENST1	ENSP1	MANE Select	NC_000001.11	100	1000	+
```

## 1) `dnakit.annotation.transcript_store_build()`

- **Purpose:** Build a gtars transcript store file from transcript records, cdot JSON files and a MANE summary.
- **Method:** `TxStoreBuilder.add_transcript` for each record (chromosome = refget digest, strand `+` or `-`), `add_chrom_mapping` for contig names, `load_cdot` for cdot files (strand 1 or -1; gzip is detected), `load_mane_summary` for an NCBI MANE summary TSV, then `build(path)`. The output file must not exist.
- **API:** `dnakit.annotation.transcript_store_build(output, *, transcripts=None, cdot_files=None, mane_summary=None, chrom_mappings=None, verbose=True)`.
- **Input:** `output`, `transcripts`, `cdot_files`, `mane_summary`, `chrom_mappings`.
- **Output:** `transcript_count`, `transcripts_added`, `cdot_transcripts_loaded`, `mane_rows_loaded`, `output`.
- **Example code:**

```python
from dnakit.annotation.gtars_p24_reftx_functions import transcript_store_build

import os
CH = 'ixSNGKM4iP6C8F3mkCr_2N0I8RB4mD5U'
def tx(acc, gene, strand, mane):
    return {'accession': acc, 'gene': gene, 'chrom': CH, 'strand': strand, 'cds_start': 150, 'cds_end': 950,
            'exons': [{'start': 100, 'end': 400}, {'start': 600, 'end': 1000}], 'mane': mane}
if os.path.exists('tx.reftx'):
    os.remove('tx.reftx')
r = transcript_store_build('tx.reftx', transcripts=[tx('NM_1.1', 'G1', '+', {'select': True, 'plus_clinical': False}),
                                                   tx('NM_2.1', 'G2', '-', None)], verbose=False)
print(r['transcripts_added'], r['transcript_count'])
if os.path.exists('cdot.reftx'):
    os.remove('cdot.reftx')
c = transcript_store_build('cdot.reftx', cdot_files=['c.json'], mane_summary='mane.tsv', chrom_mappings={'NC_000001.11': CH}, verbose=False)
print(c['cdot_transcripts_loaded'], c['mane_rows_loaded'], c['transcript_count'])
```

- **Example output:**

```text
2 2
1 2 1
```

## 2) `dnakit.annotation.transcript_record_lookup()`

- **Purpose:** Look up a transcript (gene, chromosome, strand, CDS span, exons, MANE flags) in a transcript store.
- **Method:** `ReadonlyTxStore.open`, then `lookup`, `get_chrom_accession` and `get_strand`. `found` is false when the accession is not in the store.
- **API:** `dnakit.annotation.transcript_record_lookup(transcript_store, accession, *, verbose=True)`.
- **Input:** Store path, accession (with version).
- **Output:** `found`, `transcript`, `chromosome_accession`, `strand`.
- **Example code:**

```python
from dnakit.annotation.gtars_p24_reftx_functions import transcript_record_lookup

r = transcript_record_lookup('tx.reftx', 'NM_2.1', verbose=False)
print(r['found'], r['strand'], r['chromosome_accession'])
print(r['transcript']['exons'], r['transcript']['mane'])
print(transcript_record_lookup('tx.reftx', 'NM_404.1', verbose=False)['found'])
```

- **Example output:**

```text
True - SQ.ixSNGKM4iP6C8F3mkCr_2N0I8RB4mD5U
[{'start': 100, 'end': 400}, {'start': 600, 'end': 1000}] None
False
```

## 3) `dnakit.annotation.transcript_coordinate_to_genomic()`

- **Purpose:** Convert a transcript position to a genomic position: a coding (c.) position, a non-coding (n.) position, or a c. position of a gene's MANE Select transcript.
- **Method:** `CoordinateMapper.c_to_g_full` (`coordinate_type='c'`), `n_to_g_full` (`'n'`) or `c_to_g_by_gene` (when `gene` is given; needs MANE flags in the store). `datum` is gtars' optional datum argument for c. positions.
- **API:** `dnakit.annotation.transcript_coordinate_to_genomic(transcript_store, position, *, accession=None, gene=None, coordinate_type='c', datum=None, verbose=True)`.
- **Input:** Store path, `position`, `accession` or `gene`, `coordinate_type`, `datum`.
- **Output:** `genomic_pos`, `chrom_accession`, `strand`.
- **Example code:**

```python
from dnakit.annotation.gtars_p24_reftx_functions import transcript_coordinate_to_genomic

a = transcript_coordinate_to_genomic('tx.reftx', 10, accession='NM_1.1', verbose=False)
print(a['genomic_pos'], a['strand'])
b = transcript_coordinate_to_genomic('tx.reftx', 5, accession='NM_2.1', coordinate_type='n', verbose=False)
print(b['genomic_pos'], b['strand'])
g = transcript_coordinate_to_genomic('tx.reftx', 10, gene='G1', verbose=False)
print(g['genomic_pos'], g['accession'])
```

- **Example output:**

```text
159 +
995 -
159 NM_1.1
```
