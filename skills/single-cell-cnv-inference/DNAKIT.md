# DNAKit execution: single-cell-cnv-inference

The upstream R examples use inferCNV/copyKAT/Numbat/SCEVAN and need matching counts, genomic feature coordinates and appropriate references. The output is RNA-derived evidence about large-scale DNA copy-number states, not direct DNA measurement or a validated clinical call. Retain normal-reference selection, genomic build and uncertainty in reports.

Install this complete folder with `dnakit skills install single-cell-cnv-inference --destination ./analysis-skills`. Inspect optional dependencies with `dnakit skills doctor single-cell-cnv-inference`. Resolve script/resource paths against the installed folder.

Import validation covers metadata, file integrity and Python syntax. It does not establish scientific accuracy, service availability or a successful pipeline run. Cite the upstream authors and the actual data/methods used; see UPSTREAM.md.
