# DNAKit execution: nextflow-development

The retained scripts detect input types and generate samplesheets. WGS/WES uses nf-core/sarek; ATAC uses nf-core/atacseq. Install Nextflow, its Java runtime, a supported container runtime and pipeline references separately. SRA/GEO downloads require network and SRA tools. Run scripts/check_environment.py and inspect each script --help. Pipeline test profiles and scientific workflows have not been run.

Install this complete folder with `dnakit skills install nextflow-development --destination ./analysis-skills`. Inspect optional dependencies with `dnakit skills doctor nextflow-development`. Resolve script/resource paths against the installed folder.

Import validation covers metadata, file integrity and Python syntax. It does not establish scientific accuracy, service availability or a successful pipeline run. Cite the upstream authors and the actual data/methods used; see UPSTREAM.md.
