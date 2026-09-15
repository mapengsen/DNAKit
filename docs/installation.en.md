# Installation

## Install from PyPI

Install `0.1.4`:

```bash
pip install dnakit==0.1.4
```

Package page: [pypi.org/project/dnakit](https://pypi.org/project/dnakit/)

## Other platforms

The project repository has provided the following platform adaptation files:

- [Bioconda Formula](https://github.com/mapengsen/DNAKit/tree/main/packaging/bioconda);
- [GNU Guix package definition](https://github.com/mapengsen/DNAKit/tree/main/packaging/guix);
- [Galaxy Tool Shed wrapper](https://github.com/mapengsen/DNAKit/tree/main/galaxy/dnakit).

These adaptation files are pinned to `0.1.1`, and the official installation command needs to be reviewed and released by the corresponding platform before it can be used.

## Development environment

Create and activate the complete Conda development environment in the project root directory:

```bash
conda env create -f environment-dev.yml
conda activate dnakit-dev
```

The development environment contains testing, type checking, documentation, and build tools. For additional dependencies on optional functions see
[`pyproject.toml`](https://github.com/mapengsen/DNAKit/blob/main/pyproject.toml).

## Optional bioinformatics backend

Install the optional backend used by the VEP, ClinVar, dbSNP, gnomAD, dN/dS, and Golden Gate
functions:

```bash
python -m pip install "dnakit[external-tools]"
```

The default `pip install dnakit` neither installs nor loads this backend.

## Deep-learning property prediction

Install the Enformer regulatory-track and NT Revised/Genomic Benchmarks task-checkpoint
inference dependencies:

```bash
python -m pip install "dnakit[neural,neural-enformer]"
```

The 27 task checkpoints are not distributed in the package. Prediction automatically downloads the
selected task checkpoint if absent and reuses existing files. The default cache is
`./ckpt/enformer-benchmarks/`; override it with
`checkpoint_dir` or `checkpoint_path`. See [checkpoint cache rules](api/features/23_deep_learning_property_prediction.md#checkpoint-cache); the same page describes dependencies and limitations for other models.

## Agents and MCP

Install the MCP server used by Codex, Claude, Cursor, and other MCP-capable Agents:

```bash
python -m pip install "dnakit[agent]"
```

Run `dnakit-mcp` after installation. See [LLM tools → MCP](llm_tools/mcp.md) for configuration.

# Repository

[github.com/mapengsen/DNAKit](https://github.com/mapengsen/DNAKit)

## Optional dependencies for bioinformatics Python APIs

In the current workspace, use `python -m pip install -e ".[bio]"` to install all Python analysis dependencies, or install only `.[bio-analysis]`. Executables such as MAFFT, BLAST, HMMER, BWA, samtools, and bcftools must be installed separately. See [Dependencies and implementation](api/features/27_bioinformatics.md#dependencies).
