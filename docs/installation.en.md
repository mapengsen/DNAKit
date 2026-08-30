# Install

DNAKit requires Python 3.10 or higher.

## Install from PyPI

Install `0.1.2`:

```bash
pip install dnakit==0.1.2
```

Verify installation:

```bash
python -c "import dnakit; print(dnakit.__version__)"
```

## Other platforms

The project repository has provided the following platform adaptation files:

- [Bioconda Formula](https://github.com/mapengsen/DNAKit/tree/main/packaging/bioconda);
- [GNU Guix package definition](https://github.com/mapengsen/DNAKit/tree/main/packaging/guix);
- [Galaxy Tool Shed wrapper](https://github.com/mapengsen/DNAKit/tree/main/galaxy/dnakit).

These adaptation files are fixedly corresponding to `0.1.1`, and the official installation command needs to be reviewed and released by the corresponding platform before it can be used.

## Development environment

Create and activate the complete Conda development environment in the project root directory:

```bash
conda env create -f environment-dev.yml
conda activate dnakit-dev
```

The development environment contains testing, type checking, documentation, and build tools. For additional dependencies on optional functions see
[`pyproject.toml`](https://github.com/mapengsen/DNAKit/blob/main/pyproject.toml).

# Warehouse address

[github.com/mapengsen/DNAKit](https://github.com/mapengsen/DNAKit)
