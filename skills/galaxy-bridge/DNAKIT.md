# DNAKit execution: galaxy-bridge

From this installed folder run `python galaxy_bridge.py --help`. Offline --search/--recommend/--workflow read the complete bundled catalog and curated tool resources. Actual --run requires a Galaxy account, BioBlend and the user-selected server; tool versions/availability may differ from the snapshot. The bundled demo output is an upstream example, not evidence of a newly executed analysis. Root clawbio.py commands require a full ClawBio checkout.

Install this complete folder with `dnakit skills install galaxy-bridge --destination ./analysis-skills`. Inspect optional dependencies with `dnakit skills doctor galaxy-bridge`. Resolve script/resource paths against the installed folder.

Import validation covers metadata, file integrity and Python syntax. It does not establish scientific accuracy, service availability or a successful pipeline run. Cite the upstream authors and the actual data/methods used; see UPSTREAM.md.
