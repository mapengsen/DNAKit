# DNAKit execution: bioconductor-bridge

From this installed folder run `python bioconductor_bridge.py --help`. Domain/workflow templates are available locally; live package discovery/documentation and installation require R/BiocManager/network as indicated. The small original clawbio/common/checksums.py helper is bundled under resources and the script root is adjusted to that portable location. Root `clawbio.py` commands in retained upstream prose require a separate full ClawBio checkout; use the standalone script here.

Install this complete folder with `dnakit skills install bioconductor-bridge --destination ./analysis-skills`. Inspect optional dependencies with `dnakit skills doctor bioconductor-bridge`. Resolve script/resource paths against the installed folder.

Import validation covers metadata, file integrity and Python syntax. It does not establish scientific accuracy, service availability or a successful pipeline run. Cite the upstream authors and the actual data/methods used; see UPSTREAM.md.
