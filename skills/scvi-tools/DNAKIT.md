# DNAKit execution: scvi-tools

DNA-relevant entrypoints are PeakVI for scATAC-seq and MultiVI for paired/unpaired RNA+ATAC. The complete upstream skill also contains RNA workflows; those are not counted as separate DNA additions. Install scvi-tools/PyTorch and provide count matrices, feature annotations and a suitable training environment. Training scripts are preserved; no models were trained or accuracy validated.

Install this complete folder with `dnakit skills install scvi-tools --destination ./analysis-skills`. Inspect optional dependencies with `dnakit skills doctor scvi-tools`. Resolve script/resource paths against the installed folder.

Import validation covers metadata, file integrity and Python syntax. It does not establish scientific accuracy, service availability or a successful pipeline run. Cite the upstream authors and the actual data/methods used; see UPSTREAM.md.
