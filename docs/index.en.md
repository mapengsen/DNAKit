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

DNAKit is a reproducible analysis toolkit that makes DNA and gene operations convenient. It provides:

1. **Data preparation:** DNA sequence and annotation file I/O, database downloads and queries, and validation.
2. **Data processing:** sequence standardization and editing, transcription and translation, deduplication and clustering, and dataset splitting.
3. **Data analysis:** DNA descriptors, representations and fingerprints, sequence and pattern searches, physicochemical and thermodynamic properties, and deep-learning property prediction.
4. **Data evaluation:** validity, uniqueness, diversity and novelty, plus similarity, distance, and sequence alignment.
5. **Visualization:** sequence text plots, position highlighting, linear and circular DNA maps, and alignment plots.
6. **Large language model tools:** provide a [Python API](llm_tools/python_api.md), [MCP integration](llm_tools/mcp.md), [Skills workflows](llm_tools/skills.md), and a built-in [Agent](llm_tools/agent.md) for general-purpose LLM calls and natural-language task execution.

# Support and community

**GitHub repository:** [github.com/mapengsen/DNAKit](https://github.com/mapengsen/DNAKit)

<div id="wechat-community" class="wechat-community" aria-labelledby="wechat-title">
  <strong id="wechat-title">WeChat community:</strong>
  <img src="../assets/images/community/wechat-community-qr.png" alt="QR code for the DNAKit WeChat community contact" loading="lazy">
</div>

# Cite DNAKit {#citing-dnakit}

DNAKit does not yet have a published paper or DOI. The current citation information is
maintained in the repository `README.md` and can be cited as:

```text
DNAKit contributors. DNAKit 0.1.5, 2026.
```

Citing DNAKit does not replace citations for the algorithms, parameter sets, backends, and
reference databases used in an analysis. The relevant methods, papers, databases, and websites
are listed on the [acknowledgements and primary sources page](acknowledgements.md#acknowledgements).
`Provenance`, `BackendInfo`, `ReferenceLibrary`, and `RunManifest` preserve this information in
the results.
