<p align="right">
  <a href="https://github.com/mapengsen/DNAKit/blob/main/README_CN.md">简体中文</a> | <strong>English</strong>
</p>

<p align="center">
  <img alt="DNAKit logo" src="https://raw.githubusercontent.com/mapengsen/DNAKit/main/docs/assets/images/DNAKit-icon.png" width="40%">
</p>

<p align="center">
  <a href="https://github.com/mapengsen/DNAKit/actions/workflows/ci.yml">
    <img src="https://github.com/mapengsen/DNAKit/actions/workflows/ci.yml/badge.svg?branch=main&amp;event=push" alt="CI status">
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
</p>

<p align="center">
  <strong>DNAKit: A Comprehensive Toolkit for Efficient DNA Research</strong>
</p>

# What is DNAKit?

DNAKit is a reproducible Python toolkit for DNA sequence analysis, positioned as an "RDKit for
DNA." It covers standardization, sequence and annotation formats, basic operations, descriptors,
pattern scanning, thermodynamics, fingerprints, similarity analysis, clustering and dataset
splitting, comprehensive evaluation, molecular-biology simulation, and visualization. DNA
foundation models can optionally be used to extract representations for k-means clustering.

# Installation and quick start

Install the current release:

```bash
pip install dnakit==0.1.2
```

To contribute to development, clone the repository and create the complete Conda environment:

```bash
conda env create -f environment-dev.yml
conda activate dnakit-dev
```

Then explore the [feature tree](https://mapengsen.github.io/DNAKit/en/api/features/function_tree/).

# Additional distribution platforms

The repository provides [Bioconda and GNU Guix recipes](https://github.com/mapengsen/DNAKit/tree/main/packaging)
and a [Galaxy Tool Shed wrapper](https://github.com/mapengsen/DNAKit/tree/main/galaxy/dnakit).
These files target `0.1.1`. Until they are reviewed and published by their respective platforms,
use the PyPI installation method above.

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

# Changelog

See [CHANGELOG.md](https://github.com/mapengsen/DNAKit/blob/main/CHANGELOG.md) for release notes.

# Citation

The DNAKit project paper has not yet been published. Cite the current version as:

```text
DNAKit contributors. DNAKit 0.1.2, 2026. https://github.com/mapengsen/DNAKit
```

This section will be updated after the paper is published. When using specific algorithms,
databases, or optional backends, cite their corresponding sources as well.
