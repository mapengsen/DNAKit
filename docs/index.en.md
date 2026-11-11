# DNAKit

**DNAKit: A Comprehensive Toolkit for Efficient DNA Research**

<p align="center">
  <a href="https://pypi.org/project/dnakit/">
    <img src="https://img.shields.io/pypi/v/dnakit?include_prereleases=true&amp;label=PyPI&amp;logo=pypi&amp;cacheSeconds=300" alt="PyPI version">
  </a>
  <a href="https://github.com/mapengsen/DNAKit/tree/main/packaging/bioconda">
    <img src="https://img.shields.io/badge/Bioconda-recipe-43B02A?logo=anaconda&amp;logoColor=white" alt="Bioconda recipe">
  </a>
  <a href="https://github.com/mapengsen/DNAKit/tree/main/galaxy/dnakit">
    <img src="https://img.shields.io/badge/Galaxy-wrapper-2C3143?logo=galaxy&amp;logoColor=white" alt="Galaxy wrapper">
  </a>
  <a href="https://github.com/mapengsen/DNAKit/tree/main/packaging/guix">
    <img src="https://img.shields.io/badge/GNU%20Guix-package-A42E2B?logo=gnu&amp;logoColor=white" alt="GNU Guix package">
  </a>
  <a href="llm_tools/mcp/">
    <img src="https://img.shields.io/badge/MCP-server-1A1A1A?logo=modelcontextprotocol&amp;logoColor=white" alt="MCP server">
  </a>
  <a href="llm_tools/skills/">
    <img src="https://img.shields.io/badge/Skills-agent%20workflows-D97757?logo=data:image/svg%2bxml;base64,PHN2ZyB4bWxucz0iaHR0cDovL3d3dy53My5vcmcvMjAwMC9zdmciIHZpZXdCb3g9IjAgMCAyNCAyNCI+PHBhdGggZmlsbD0id2hpdGUiIGQ9Ik0xMyAxIDMgMTRoN2wtMS41IDlMMjAgOS41aC03LjJ6Ii8+PC9zdmc+" alt="Agent Skills">
  </a>
  <a href="https://dnakit.mapengsen.asia">
    <img src="https://img.shields.io/badge/Web%20App-online-2EA44F?logo=googlechrome&amp;logoColor=white" alt="Web app">
  </a>
  <a href="gui/">
    <img src="https://img.shields.io/badge/GUI-Windows-0078D6?logo=data:image/svg%2bxml;base64,PHN2ZyB4bWxucz0iaHR0cDovL3d3dy53My5vcmcvMjAwMC9zdmciIHZpZXdCb3g9IjAgMCAyNCAyNCI+PHBhdGggZmlsbD0id2hpdGUiIGQ9Ik0wIDBoMTEuNHYxMS40SDB6TTEyLjYgMEgyNHYxMS40SDEyLjZ6TTAgMTIuNmgxMS40VjI0SDB6TTEyLjYgMTIuNkgyNFYyNEgxMi42eiIvPjwvc3ZnPg==" alt="GUI for Windows">
  </a>
  <a href="gui/">
    <img src="https://img.shields.io/badge/GUI-macOS-000000?logo=apple&amp;logoColor=white" alt="GUI for macOS">
  </a>
  <a href="gui/">
    <img src="https://img.shields.io/badge/GUI-Linux-FCC624?logo=linux&amp;logoColor=black" alt="GUI for Linux">
  </a>
</p>

DNAKit is a unified toolkit for DNA analysis. It brings the functions of 284 existing software packages, web services and agent tool collections under one set of conventions and provides 6,698 operations, so researchers and AI agents can use the same functions for a whole analysis, from data download to model training.

<div class="dnakit-stats">
  <div><strong>6,698</strong><span>operations</span></div>
  <div><strong>12</strong><span>DNA analysis domains</span></div>
  <div><strong>148</strong><span>file formats</span></div>
  <div><strong>107</strong><span>public databases</span></div>
  <div><strong>48</strong><span>DNA foundation models</span></div>
  <div><strong>7</strong><span>interfaces</span></div>
</div>

# Why DNAKit

- **One set of conventions, no per-tool adapters:** every operation uses the same naming, data formats, return values and error messages, so you no longer write glue code between tools.
- **Install once, use it seven ways:** the [Python API](llm_tools/python_api.md), the command line, the [web app](https://dnakit.mapengsen.asia), the [desktop app](gui.md), [MCP](llm_tools/mcp.md), [Skills](llm_tools/skills.md) and the built-in [Agent](llm_tools/agent.md) share one implementation and give the same results.
- **DNA tools ready for AI agents:** 6,120 operations can be called through MCP, and 575 workflow Skills are included; the same request gives the same answer in Codex, Claude Code and Gemini CLI.
- **No coding required:** the web app runs in a browser with nothing to install; the desktop app (`dnakit gui`) runs on Windows, macOS and Linux.
- **Every step is traceable:** each operation is one explicit function call, and its inputs, parameters, backend and reference database are recorded so that an analysis can be checked and reproduced.
- **Extend it yourself:** a Python function with type annotations and a docstring is registered as a tool; you can also describe what you need in plain language and let the Agent generate the tool, which is registered after it passes the checks.

# What it can do

| Domain | Operations | Examples |
|:---|---:|:---|
| Data & databases | 692 | Reading, writing and converting 148 file formats such as FASTA, BAM, VCF and BED; querying and downloading from 107 databases such as NCBI, Ensembl, ENCODE and ClinVar |
| Sequence analysis | 251 | Reverse complement, transcription and translation, ORF finding, k-mers and sequence composition, DNA descriptors |
| Alignment & homology | 328 | Sequence alignment, similarity and distance, BLAST homology search |
| Reads & assembly | 566 | Read trimming and quality control, UMI handling, read mapping, genome assembly |
| Variants & clinical | 907 | Variant calling, VEP annotation, ClinVar clinical significance, structural variants |
| Population & GWAS | 680 | Population-genetic statistics, haplotypes and linkage disequilibrium, genome-wide association analysis |
| Genes & expression | 1,009 | Gene annotation and orthologs, expression and enrichment analysis, single-cell and spatial transcriptomics |
| Regulation & epigenome | 666 | Motif scanning, peak analysis, methylation, Hi-C contact matrices |
| Evolution & microbiome | 452 | Phylogenetic trees, dN/dS, pangenomes, metagenomics and taxonomy |
| Design & engineering | 389 | Primer design, melting temperature and thermodynamics, restriction digestion, Golden Gate assembly, CRISPR design |
| Sequence machine learning | 195 | Dataset deduplication, clustering and splitting, fingerprints and features, training property-prediction models |
| DNA foundation models | 233 | Representation extraction, fine-tuning and prediction with 48 models, including AlphaGenome and Evo 2 |

DNAKit also evaluates generated sequences (validity, uniqueness, diversity and novelty) and draws figures (sequence text plots, position highlighting, linear and circular DNA maps, and alignment plots). A further 330 generic operations, such as infrastructure and literature search, fall outside these domains. The full list is in the [function tree](api/features/function_tree.md).

# Quick start

```bash
pip install dnakit
```

```python
from dnakit.ops import reverse_complement

print(reverse_complement("AACG"))  # CGTT
```

```bash
pip install "dnakit[agent]"
# register the MCP server
claude mcp add dnakit -- dnakit-mcp     # Claude Code
codex mcp add dnakit -- dnakit-mcp      # Codex
# install the Skills (or --agent claude-code, --agent cursor)
dnakit skills install dnakit --agent codex
```

See [Installation](installation.md) for other ways to install, and [MCP](llm_tools/mcp.md) and [Skills](llm_tools/skills.md) for agent use.

# Support and community

**GitHub repository:** [github.com/mapengsen/DNAKit](https://github.com/mapengsen/DNAKit)

# Cite DNAKit {#citing-dnakit}

DNAKit does not yet have a published paper or DOI. The current citation information is
maintained in the repository `README.md` and can be cited as:

```text
DNAKit contributors. DNAKit 0.2.0, 2026.
```

Citing DNAKit does not replace citations for the algorithms, parameter sets, backends, and
reference databases used in an analysis. The relevant methods, papers, databases, and websites
are listed on the [acknowledgements and primary sources page](acknowledgements.md#acknowledgements).
`Provenance`, `BackendInfo`, `ReferenceLibrary`, and `RunManifest` preserve this information in
the results.
