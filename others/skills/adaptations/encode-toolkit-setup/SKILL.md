---
name: encode-toolkit-setup
description: Set up ENCODE access through DNAKit, including MCP discovery, public queries, credentials and connection troubleshooting.
disable-model-invocation: true
---

## DNAKit execution

Use `dnakit skills doctor encode-toolkit-setup` to inspect scientific dependencies. Read [DNAKIT.md](DNAKIT.md) for this skill's DNAKit functions and bundled script entrypoints. Keep the workflow, input checks and reporting requirements below.


<!-- DNAKit metadata portability edit, 2026-09-12; see UPSTREAM.md. -->

# ENCODE Toolkit Setup

## When to Use

- User needs help installing or configuring the ENCODE Toolkit MCP server
- User is getting connection errors or server startup failures
- User asks "how do I set up ENCODE?" or "install ENCODE toolkit"
- User needs to configure ENCODE credentials for restricted data access
- User wants to verify their ENCODE server connection is working
- User is setting up a new environment and needs the ENCODE plugin

Help the user set up the ENCODE Toolkit server. The server connects Claude to the ENCODE Project genomics database — the largest public catalog of functional genomic elements with 8,000+ experiments across 50+ assay types.

## Installation and connection

Install DNAKit in the same Python environment that the Agent will launch:

```bash
python -m pip install "dnakit[agent]"
dnakit skills install encode-toolkit-setup --agent codex
codex mcp add dnakit -- /absolute/path/to/dnakit-mcp
codex mcp list
```

Use the actual executable path from that environment. Restart the client after
changing its MCP configuration. Clients using JSON `mcpServers` configuration
can register the same executable as `{"mcpServers":{"dnakit":{"command":"/absolute/path/to/dnakit-mcp"}}}`.

DNAKit exposes compact MCP discovery tools. Search with `search_dnakit_tools`,
inspect the chosen function with `describe_dnakit_tool`, then execute it with
`call_dnakit_tool`. ENCODE functions are named `dnakit.encode.<operation>`.

Verify a public query before downloads or tracking:

```python
import json
import dnakit

result = json.loads(dnakit.encode.search_experiments(
    assay_title="ATAC-seq", organism="Homo sapiens", limit=3,
))
print(result)
```

Check returned errors and actual experiment accessions; do not treat a registered
MCP configuration as proof that a network query succeeded. Verify filter values
with `dnakit.encode.get_metadata` and data availability with `dnakit.encode.get_facets`.

## Authentication

Install `dnakit[encode-credentials]` when using persistent credential storage. DNAKit stores its ENCODE tracker and fallback credentials under `~/.dnakit-encode/`.

Most ENCODE data is public and needs no authentication. For restricted/unreleased data:

1. Get API credentials from https://www.encodeproject.org/profile/ (requires ENCODE account)
2. Store them:
   ```
   Ask: "Store my ENCODE credentials"
   → Calls dnakit.encode.manage_credentials(action="store", access_key="...", secret_key="...")
   ```
3. Credentials are encrypted via the OS keyring (macOS Keychain, Linux Secret Service)
4. To verify: `dnakit.encode.manage_credentials(action="status")`
5. To remove: `dnakit.encode.manage_credentials(action="remove")`

---

## 20 Available Tools

After setup, these tools are available:

| Category | Tools | Purpose |
|----------|-------|---------|
| **Search** | `dnakit.encode.search_experiments`, `dnakit.encode.get_facets`, `dnakit.encode.get_metadata` | Find experiments, explore data landscape, get valid filter values |
| **Experiment Details** | `dnakit.encode.get_experiment`, `dnakit.encode.compare_experiments` | Get full experiment metadata, compare two experiments |
| **Files** | `dnakit.encode.search_files`, `dnakit.encode.list_files`, `dnakit.encode.get_file_info` | Find files, list files for an experiment, get file details |
| **Download** | `dnakit.encode.download_files`, `dnakit.encode.batch_download` | Download individual or batch files with MD5 verification |
| **Tracking** | `dnakit.encode.track_experiment`, `dnakit.encode.list_tracked`, `dnakit.encode.summarize_collection` | Local experiment tracking with SQLite |
| **Provenance** | `dnakit.encode.log_derived_file`, `dnakit.encode.get_provenance` | Log analysis outputs with full lineage |
| **Citations** | `dnakit.encode.get_citations`, `dnakit.encode.link_reference` | Publication data, cross-reference to PubMed/GEO |
| **Credentials** | `dnakit.encode.manage_credentials` | Store/remove API credentials |
| **Collection** | `dnakit.encode.summarize_collection` | Summarize tracked experiment portfolio |

---

## First-Run Walkthrough: Pancreatic Islet Epigenomics

This walkthrough demonstrates a complete workflow from installation to data exploration.

### 1. Explore what's available
```
"What ENCODE assay types are available for human pancreas?"
→ dnakit.encode.get_facets(facet_field="assay_title", organ="pancreas", organism="Homo sapiens")
```

### 2. Find specific experiments
```
"Find all histone ChIP-seq experiments on human pancreas"
→ dnakit.encode.search_experiments(assay_title="Histone ChIP-seq", organ="pancreas", organism="Homo sapiens")
```

### 3. Examine an experiment
```
"Get details for ENCSR123ABC"
→ dnakit.encode.get_experiment(accession="ENCSR123ABC")
```

### 4. Find the right files
```
"List the preferred BED files for ENCSR123ABC"
→ dnakit.encode.list_files(accession="ENCSR123ABC", file_format="bed", assembly="GRCh38")
```

### 5. Download data
```
"Download the IDR-thresholded peaks for ENCSR123ABC"
→ dnakit.encode.download_files(accession="ENCSR123ABC", file_format="bed", output_type="IDR thresholded peaks")
```

### 6. Track your experiment
```
"Track ENCSR123ABC in my local database with note 'H3K27ac pancreatic islets'"
→ dnakit.encode.track_experiment(accession="ENCSR123ABC", notes="H3K27ac pancreatic islets")
```

---

## Cross-Database Integration

The ENCODE Toolkit works alongside other MCP servers and REST APIs:

| Database | Access Method | What It Adds |
|----------|--------------|--------------|
| **PubMed** | MCP server (`search_articles`) | Literature citations for ENCODE experiments |
| **bioRxiv** | MCP server (`search_preprints`) | Preprint discovery for latest research |
| **ClinicalTrials.gov** | MCP server (`search_trials`) | Clinical trial cross-reference |
| **Open Targets** | MCP server (`query_open_targets_graphql`) | Drug target identification |
| **GTEx** | REST API via skill | Tissue-specific expression context |
| **ClinVar** | REST API via skill | Clinical variant annotation |
| **GWAS Catalog** | REST API via skill | Trait-associated variant lookups |
| **gnomAD** | GraphQL via skill | Population allele frequencies |
| **Ensembl** | REST API via skill | VEP annotation, Regulatory Build |
| **UCSC** | REST API via skill | Genome browser tracks, cCRE data |
| **GEO** | E-utilities via skill | Complementary expression datasets |
| **JASPAR** | REST API via skill | TF binding motif databases |
| **CellxGene** | REST API via skill | Single-cell expression atlases |

---

## 47 Expert Skills

Beyond the 20 tools, the ENCODE Toolkit includes 47 skills providing domain expertise:

- **Core (5)**: setup, search-encode, download-encode, track-experiments, cross-reference
- **Analysis (9)**: quality-assessment, integrative-analysis, regulatory-elements, epigenome-profiling, compare-biosamples, visualization-workflow, motif-analysis, peak-annotation, batch-analysis
- **Pipelines (7)**: pipeline-chipseq, pipeline-atacseq, pipeline-rnaseq, pipeline-wgbs, pipeline-hic, pipeline-dnaseseq, pipeline-cutandrun
- **External DBs (9)**: gtex-expression, clinvar-annotation, cellxgene-context, gwas-catalog, jaspar-motifs, ensembl-annotation, geo-connector, gnomad-variants, ucsc-browser
- **Workflows (10)**: data-provenance, cite-encode, variant-annotation, pipeline-guide, single-cell-encode, disease-research, publication-trust, bioinformatics-installer, scientific-writing, liftover-coordinates
- **Data Aggregation (4)**: histone-aggregation, accessibility-aggregation, hic-aggregation, methylation-aggregation
- **Meta-Analysis (2)**: scrna-meta-analysis, multi-omics-integration
- **Functional Genomics (1)**: functional-screen-analysis

---

## Pitfalls & Troubleshooting

| Problem | Cause | Fix |
|---------|-------|-----|
| "Server not found" | Claude not restarted after config change | Restart Claude Desktop / reload Claude Code |
| "uvx not found" | uv not installed | `curl -LsSf https://astral.sh/uv/install.sh \| sh` |
| Timeout errors | Slow connection or ENCODE API load | Retry; rate limit (10 req/sec) is handled automatically |
| 403 on downloads | File requires authentication | `dnakit.encode.manage_credentials(action="store", ...)` |
| No results returned | Filters too narrow | Broaden filters; use `dnakit.encode.get_facets` to see available data |
| "Invalid accession" | Wrong format | Must be ENCSR/ENCFF/ENCBS format (e.g., ENCSR000AAA) |
| Empty facets | API connectivity issue | Check internet; try `dnakit.encode.get_metadata(metadata_type="assays")` |
| Stale results | Cached data | Cache TTL is 1 hour; restart server to clear |

---

## Code Examples

### 1. Verify server connection with metadata query
```
dnakit.encode.get_metadata(metadata_type="assays")
```

Expected output:
```json
{
  "assays": ["ATAC-seq", "ChIP-seq", "CUT&RUN", "CUT&Tag", "DNase-seq", "Hi-C", "MPRA", "RNA-seq", "STARR-seq", "WGBS", "eCLIP", "scATAC-seq", "scRNA-seq"]
}
```

### 2. Test search functionality
```
dnakit.encode.search_experiments(assay_title="ATAC-seq", organ="brain", organism="Homo sapiens", limit=3)
```

Expected output:
```json
{
  "total": 32,
  "results": [
    {"accession": "ENCSR000AAA", "assay_title": "ATAC-seq", "biosample_summary": "brain", "status": "released"}
  ]
}
```

### 3. Test facet exploration
```
dnakit.encode.get_facets(facet_field="organ", organism="Homo sapiens")
```

Expected output:
```json
{
  "facets": {
    "organ": {"brain": 450, "blood": 380, "liver": 220, "heart": 180, "lung": 150}
  }
}
```

## Related Skills

| Skill | When to Use |
|-------|------------|
| `search-encode` | First skill to use after setup — find experiments by assay, tissue, target |
| `download-encode` | Download ENCODE files (BED, bigWig, FASTQ, BAM) after finding experiments |
| `pipeline-guide` | Set up Nextflow pipelines for processing raw ENCODE data |
| `bioinformatics-installer` | Install all bioinformatics tools needed for ENCODE analysis |
| `cross-reference` | Link ENCODE experiments to PubMed, GEO, ClinicalTrials.gov |
| `quality-assessment` | Evaluate data quality before analysis |
| `publication-trust` | Verify literature claims backing analytical decisions |

---

## Presenting Results

When reporting setup results:

- **Connection status**: Confirm the ENCODE Toolkit server is connected and responding. Report the server version if available
- **Available tools**: List the 20 available ENCODE tools grouped by function (search, download, track, cross-reference, credentials)
- **Test query result**: Run a simple validation query (e.g., `dnakit.encode.get_metadata(metadata_type="assays")`) and confirm it returns results successfully
- **Authentication status**: Note whether credentials are configured (for restricted data) or that public data access requires no authentication
- **Troubleshooting**: If any issues were encountered during setup, summarize the problem and resolution
- **Next steps**: Suggest `search-encode` to find experiments, or `dnakit.encode.get_facets` to explore what ENCODE data is available for their research area

## For the request: "$ARGUMENTS"
