# Calling conventions and examples for the added APIs

The [complete index](index.md) links all **484 added public functions**, with individual signatures, parameter descriptions and returns. These APIs work directly in Python without installing a Skill.

## Install the source checkout

These changes have not yet been published to PyPI. From the repository root:

```bash
python -m pip install -e "."
```

No ToolUniverse initialization is required. Scientific programs, model weights, reference data and API credentials remain task-specific requirements. See [third-party notices](../../acknowledgements.md#third-party-notices).

## Discover APIs offline

```python
import dnakit

matches = dnakit.databases.search("gnomad", limit=5, verbose=False)
print([item["name"] for item in matches])
info = dnakit.databases.describe("gnomad_get_sv_by_region", verbose=False)
print(info["input_schema"])
print(len(dnakit.databases.list_operations(verbose=False)))  # 450 in this version
```

These calls read the bundled catalog without accessing a provider. The 450 operations include the `execute_tool` compatibility dispatcher. Five additional catalog/calling helpers give 455 public database functions; this does not count independent scientific algorithms.

```bash
dnakit database search "gnomad"
dnakit database info gnomad_get_sv_by_region
```

You can also search the website for a function such as [gnomad_get_sv_by_region](variants.md#dnakit.databases.gnomad_get_sv_by_region).

## Query a database

This example accesses Ensembl; its response depends on the live service and query date.

```python
import dnakit

gene = dnakit.databases.ensembl_lookup_gene(id="BRCA1", species="homo_sapiens")
print(gene)
```

For dynamic dispatch, choose one of these equivalent forms:

```python
import dnakit

params = {"id": "BRCA1", "species": "homo_sapiens"}
result = dnakit.databases.call("ensembl_lookup_gene", params)
# result = dnakit.databases.execute({"name": "ensembl_lookup_gene", "arguments": params})
# result = dnakit.databases.execute_tool(tool_name="ensembl_lookup_gene", arguments=params)
print(result)
```

`describe()["input_schema"]` gives the actual DNAKit input constraints; `config.parameter` preserves the provider's original schema. A generated wrapper's default `None` normally omits the parameter and does not specify the provider's default. Respect required fields, aliases and conditional constraints.

Results retain provider-specific envelopes and pagination; do not assume `result["data"]` is always present. Registered result schemas are structural references, and extended query branches may return different structures. Provider errors and oversized results raise exceptions. `call(..., max_result_bytes=...)` limits serialized results to 20 MB by default.

Coordinate and assembly conventions vary by operation. ChIP-Atlas `requires_data` and `requires_submission` statuses indicate missing data or submission work, not completed analysis; see [regulation](regulation.md).

## Search and download a sequence

This uses the extended existing `entrez` plus the new `accessions` and NCBI `sequence` functions. It queries a remote service and writes a file in the current directory.

```python
import dnakit

found = dnakit.search.ncbi.entrez(
    "nuccore", organism="Escherichia coli", strain="K-12",
    seq_type="complete_genome", retmax=3,
)
accessions = dnakit.search.ncbi.accessions(found.metadata["uids"])
if accessions:
    downloaded = dnakit.download.ncbi.sequence(accessions[0], "ecoli_k12.fasta", format="fasta")
    print(downloaded.files[0].path)
    print(downloaded.metadata["sequence_length"])
    print(downloaded.manifest_path)
else:
    print("No sequence matches; review the query.")
```

For ENA use an INSDC accession, for example `dnakit.download.ena.sequence("U00096.3", "ecoli_ena.fasta")`. Download functions save files and manifests, whereas database sequence wrappers return provider records. See [formats, dependencies and coordinates](../features/15_download.md#accession-sequences).

## ENCODE queries and download previews

```python
import json
import dnakit

experiments = json.loads(dnakit.encode.search_experiments(
    assay_title="ATAC-seq", organism="Homo sapiens", limit=3,
))
print(experiments)
preview = json.loads(dnakit.encode.batch_download(
    download_dir="downloads/encode", assay_title="ATAC-seq",
    assembly="GRCh38", file_format="bed", limit=3, dry_run=True,
))
print(preview)
```

`dry_run=True` queries the file list without downloading the files. Set it to `False` to download. The [20 ENCODE functions](encode.md) return text: mostly JSON, but citation and table exports also support BibTeX, RIS, CSV and TSV. Do not parse those formats with `json.loads()`. Tracking, credentials and downloads can write local state.

## Run a bundled coordinate conversion

This runs a packaged Python script and converts BED 0-based half-open coordinates to GFF 1-based inclusive coordinates, not between genome assemblies.

```python
import dnakit

print(dnakit.workflows.scripts("genomic-coordinates", verbose=False))
print(dnakit.workflows.resource_path(
    "genomic-coordinates", "scripts/convert_coords.py", verbose=False,
))
result = dnakit.workflows.run(
    "genomic-coordinates", "scripts/convert_coords.py",
    arguments=["--from", "bed", "--to", "gff", "chr1", "999", "1000", "--format", "json"],
    timeout_seconds=60, verbose=False,
)
if result["returncode"] != 0:
    raise RuntimeError(result["output"])
print(result["output"])
```

`run_support()` accepts `clinical-variant-reporter`, `clinical-trial-finder`, `omics-target-evidence-mapper`, `gwas-lookup` and `fine-mapping`. Use `arguments=["--help"]` to inspect their CLI arguments and prepare the scientific dependencies imported by each script. See [workflow parameters and resource verification](execution.md).

## Simulate independent-locus inheritance

This synthetic locus needs no remote service. With A/A and C/C parents and no mutation, each offspring receives A/C.

```python
import dnakit

locus = {
    "chromosome": "1", "position": 1, "ref": "A", "alt": "C",
    "dominance": "additive", "effect_size": 1,
}
father = {"id": "father", "sex": "Male", "loci": {"L1": {**locus, "alleles": ["A", "A"]}}}
mother = {"id": "mother", "sex": "Female", "loci": {"L1": {**locus, "alleles": ["C", "C"]}}}
traits = {"traits": {"dosage": {"loci": [dict(locus, id="L1", effect=1)]}}}
result = dnakit.workflows.simulate_inheritance(
    father, mother, traits, {"diseases": {}},
    offspring_count=2, mutation_rate=0, seed=13, verbose=False,
)
print([child["loci"]["L1"]["alleles"] for child in result["offspring"]])
# [['A', 'C'], ['A', 'C']]
```

This is independent-locus segregation and mutation, not chromosome crossover or clinical outcome prediction. See [all parameters](analysis.md#dnakit.workflows.simulate_inheritance).

## Render an existing report

```bash
python -m pip install -e ".[reports]"
```

```python
from pathlib import Path
from tempfile import TemporaryDirectory
import dnakit

with TemporaryDirectory() as directory:
    root = Path(directory)
    report = root / "report.md"
    report.write_text("# Results\nUser-supplied report content.\n", encoding="utf-8")
    result = dnakit.workflows.write_wes_report(report, root / "pdf", verbose=False)
    print(result["output"])
    print(result["manifest"])
    print(Path(result["output"]).exists())  # True
```

Use your own output directory to keep the files. This renderer preserves supplied content without generating medical interpretation. `language` selects English or Spanish labels without translating the text. See [title, logos, font and overwrite options](execution.md#dnakit.workflows.write_wes_report).

## MCP and validation scope

The same functions are available through tool discovery, description and invocation in the [MCP guide](../../llm_tools/mcp.md). Python, CLI and MCP exposures do not count as separate capabilities.

Documentation checks verify 484 unique function entries, current signatures and parameters, and build both languages. Documentation coverage and offline example checks do not establish end-to-end validation of all remote services, models or scientific workflows.
