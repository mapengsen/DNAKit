# Installation

## Install from PyPI

Install `0.1.5`:

```bash
pip install dnakit==0.1.5
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
`checkpoint_dir` or `checkpoint_path`. See [prediction dependencies and checkpoint caching](faq/tutorial_integrations.md#property-prediction). Task examples are organized into [regulatory tracks](api/features/23_prediction_regulatory_tracks.md), [epigenetic classification](api/features/23_prediction_chromatin_classification.md), and [promoter and enhancer classification](api/features/23_prediction_promoter_enhancer.md).

## Agents and MCP

Install the MCP server used by Codex, Claude, Cursor, and other MCP-capable Agents:

```bash
python -m pip install "dnakit[agent]"
```

Run `dnakit-mcp` after installation. See [LLM tools → MCP](llm_tools/mcp.md) for configuration.

# Repository

[github.com/mapengsen/DNAKit](https://github.com/mapengsen/DNAKit)

## Optional dependencies for bioinformatics Python APIs

Install the current source checkout's common genomics capabilities with Python ≥3.10:

```bash
python -m pip install -e ".[genomics]"
```

The `genomics` extra installs NumPy, Pyrodigal, pysam and PyFastANI for genotype PCA, genomic relationships and breeding prediction, prokaryotic CDS prediction, BAM statistics and genome ANI. These APIs do not need separate PLINK, Prodigal, samtools or FastANI executables. See the [genomics tutorial](api/features/38_genomics_extensions.md#pip-genomics) for their scope and examples.

These additions are not yet published to PyPI; after publication, use `python -m pip install "dnakit[genomics]"`. Wheel availability depends on the dependency and platform; platforms without a matching wheel may need a compiler. This extra does not install all 43 external analyses, reference databases or model weights.

In the current workspace, use `python -m pip install -e ".[bio]"` for general bioinformatics Python dependencies, or install only `.[bio-analysis]`. Use `.[genomics]` above for the new genomics computing dependencies. External command adapters for MAFFT, BLAST, HMMER, BWA, samtools and bcftools still require their executables. See [Dependencies and implementation](faq/bioinformatics.md#dependencies).

## Install the source checkout {#source-python-usage}

These changes have not yet been published to PyPI. From the repository root:

```bash
python -m pip install -e "."
```

No ToolUniverse initialization is required. Scientific programs, model weights, reference data and API credentials remain task-specific requirements. See [third-party notices](acknowledgements.md#third-party-notices).

## Optional genomics extensions

From this checkout, install `.[genomics-genes]` for Pyrodigal, `.[dna-storage]` for Chamaeleo, or `.[dna-nanotechnology]` for scadnano. Breeding uses `.[bio-analysis]`; the MCP client uses `.[agent]`.

Install the local gene-finding extension:

```bash
python -m pip install -e ".[genomics-genes]"
```

External executables, references and model weights remain separate; see [external genomics and MCP](api/features/39_external_genomics.md).
