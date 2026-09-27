# DNAKit execution: galaxy-bioblend

This is a transparent SKILL.md entrypoint around the original galaxy.md REST and galaxy-bioblend.md SDK guides. Use `from bioblend.galaxy import GalaxyInstance` with the user-selected server. History/dataset/workflow primitives complement the separately packaged galaxy-bridge recommendation CLI. references/original/setup-credentials.sh is retained as an optional manual helper and was not executed. The unrelated gene-mapping-by-log-fold-change material is excluded because numerical expression similarity cannot establish gene identity.

Install this complete folder with `dnakit skills install galaxy-bioblend --destination ./analysis-skills`. Inspect optional dependencies with `dnakit skills doctor galaxy-bioblend`. Resolve script/resource paths against the installed folder.

Import validation covers metadata, file integrity and Python syntax. It does not establish scientific accuracy, service availability or a successful pipeline run. Cite the upstream authors and the actual data/methods used; see UPSTREAM.md.
