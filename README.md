<p align="right">
  <a href="https://github.com/mapengsen/DNAKit/blob/main/README_CN.md">简体中文</a> | <strong>English</strong>
</p>

<p align="center">
  <img alt="DNAKit logo" src="https://raw.githubusercontent.com/mapengsen/DNAKit/main/docs/assets/images/DNAKit-icon.png" width="40%">
</p>

<p align="center">
  <a href="https://github.com/mapengsen/DNAKit/actions/workflows/ci.yml">
    <img src="https://github.com/mapengsen/DNAKit/actions/workflows/ci.yml/badge.svg?branch=main&event=push" alt="CI status">
  </a>
  <a href="https://mapengsen.github.io/DNAKit/en/">
    <img src="https://img.shields.io/badge/docs-GitHub%20Pages-4051b5" alt="Documentation">
  </a>
  <a href="https://www.python.org/downloads/">
    <img src="https://img.shields.io/badge/Python-3.10%2B-3776AB" alt="Python 3.10+">
  </a>
  <a href="https://github.com/mapengsen/DNAKit/blob/main/LICENSE">
    <img src="https://img.shields.io/github/license/mapengsen/DNAKit" alt="MIT license">
  </a>
</p>

<p align="center">
  <a href="https://pypi.org/project/dnakit/">
    <img src="https://img.shields.io/pypi/v/dnakit?include_prereleases=true&label=PyPI&logo=pypi&cacheSeconds=300" alt="PyPI version">
  </a>
  <a href="https://github.com/mapengsen/DNAKit/tree/main/packaging/bioconda">
    <img src="https://img.shields.io/badge/Bioconda-recipe-43B02A?logo=anaconda&logoColor=white" alt="Bioconda recipe">
  </a>
  <a href="https://github.com/mapengsen/DNAKit/tree/main/galaxy/dnakit">
    <img src="https://img.shields.io/badge/Galaxy-wrapper-2C3143?logo=galaxy&logoColor=white" alt="Galaxy wrapper">
  </a>
  <a href="https://github.com/mapengsen/DNAKit/tree/main/packaging/guix">
    <img src="https://img.shields.io/badge/GNU%20Guix-package-A42E2B?logo=gnu&logoColor=white" alt="GNU Guix package">
  </a>
</p>

<p align="center">
  <strong>DNAKit: A Comprehensive Toolkit for Efficient DNA Research</strong>
</p>

# What is DNAKit?

DNAKit is a reproducible Python toolkit for DNA sequence analysis, positioned as an "RDKit for
DNA." It covers standardization, sequence and annotation formats, basic operations, descriptors,
pattern scanning, thermodynamics, fingerprints, similarity analysis, clustering and dataset
splitting, comprehensive evaluation, molecular-biology simulation, and visualization. DNA
foundation models can optionally be used for direct property prediction or to extract
representations for k-means clustering.

# Installation and quick start

Install the current release:

```bash
pip install dnakit==0.1.5
```

Install the optional VEP, ClinVar, dbSNP, gnomAD, dN/dS, and Golden Gate
backend with `pip install "dnakit[external-tools]"`.

To expose DNAKit to an MCP-capable Agent, run `pip install "dnakit[agent]"` and configure the
`dnakit-mcp` command in the Agent client.

This source checkout's built-in `run_agent()` also supports optional tool generation, description optimization or compatible implementation repairs, local file result caching/resume, and retrieval from a user-selected knowledge directory (RAG). See the [Agent tutorial](docs/llm_tools/agent.en.md#agent-actions) for configuration and examples; generation and optimization require a writable checkout and explicit permission for side effects.

This source checkout also includes a curated DNA/genomics Skill catalog. Install the checkout
with `python -m pip install -e "."`, then search with `dnakit skills search "variant annotation"`
or install the discovery skill with `dnakit skills install dnakit --agent codex`.
See the [Skills tutorial](docs/llm_tools/skills.en.md) for Claude Code, Cursor, MCP, and citations.
DNAKit code remains MIT; bundled Skills retain their [original licenses and attribution](THIRD_PARTY_NOTICES.md).

This checkout adds [breeding, gene finding, methylation, DNA storage and nanostructure APIs](docs/api/features/38_genomics_extensions.en.md), plus [external analyses and MCP connections](docs/api/features/39_external_genomics.en.md).
Install selected extras locally, for example `python -m pip install -e ".[bio-analysis,genomics-genes,agent]"`.
Original projects, authors and remaining integration gaps are recorded in the [source FAQ](docs/faq/genomics_sources.en.md). The [335-project review](docs/faq/user_catalog.en.md) records every supplied project, duplicate, added interface and deferred dependency.

To contribute to development, clone the repository and create the complete Conda environment:

```bash
conda env create -f environment-dev.yml
conda activate dnakit-dev
```

# Documentation

The complete English documentation is available from the
[DNAKit documentation homepage](https://mapengsen.github.io/DNAKit/en/), including the
[feature tree](https://mapengsen.github.io/DNAKit/en/api/features/function_tree/) and
[FAQ](https://mapengsen.github.io/DNAKit/en/faq/).

# Support and community

**GitHub repository:**

[github.com/mapengsen/DNAKit](https://github.com/mapengsen/DNAKit)

For questions, feedback, or suggestions, see:

- [FAQ](https://mapengsen.github.io/DNAKit/en/faq/)
- [GitHub Issues](https://github.com/mapengsen/DNAKit/issues)

# Additional distribution platforms

The repository provides [Bioconda and GNU Guix recipes](https://github.com/mapengsen/DNAKit/tree/main/packaging)
and a [Galaxy Tool Shed wrapper](https://github.com/mapengsen/DNAKit/tree/main/galaxy/dnakit).
These files target `0.1.1`. Until they are reviewed and published by their respective platforms,
use the PyPI installation method above.

# Changelog

See [CHANGELOG.md](https://github.com/mapengsen/DNAKit/blob/main/CHANGELOG.md) for release notes.

# Citation

The DNAKit project paper has not yet been published. Cite the current version as:

```text
DNAKit contributors. DNAKit 0.1.5, 2026. https://github.com/mapengsen/DNAKit
```

This section will be updated with the paper citation after publication. When using specific algorithms, databases, or optional backends, cite their corresponding sources as well.

The working tree also provides native DNA encodings, generative metrics, statistics and specification-based optimization, plus optional format codecs and specialist alignment/NGS backends. Usage is organized by task in the [feature tutorials](docs/api/features/function_tree.en.md); see [bioinformatics notes](docs/faq/bioinformatics.en.md) for dependencies and implementation details. Install these with `python -m pip install -e ".[bio]"`.
