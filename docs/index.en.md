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
</p>

DNAKit is a reproducible toolkit for DNA sequence analysis. It provides unified objects,
sequence and annotation I/O, standardization, sequence operations, descriptors, pattern
scanning, thermodynamics, fingerprints, similarity analysis, dataset preparation,
comprehensive evaluation, molecular-biology simulation, and visualization. DNA foundation
models can optionally be used to extract representations for k-means clustering. DNAKit does
not include task-specific deep-learning predictors for promoter activity, expression levels,
transcription-factor binding strength, or CRISPR editing efficiency.

# Support and community

**GitHub repository:** [github.com/mapengsen/DNAKit](https://github.com/mapengsen/DNAKit)

# Cite DNAKit {#citing-dnakit}

DNAKit does not yet have a published paper or DOI. The current citation information is
maintained in the repository `README.md` and can be cited as:

```text
DNAKit contributors. DNAKit 0.1.3, 2026.
```

Citing DNAKit does not replace citations for the algorithms, parameter sets, backends, and
reference databases used in an analysis. The relevant methods, papers, databases, and websites
are listed on the [acknowledgements and primary sources page](acknowledgements.md#acknowledgements).
`Provenance`, `BackendInfo`, `ReferenceLibrary`, and `RunManifest` preserve this information in
the results.
