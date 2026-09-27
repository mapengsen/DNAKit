# DNAKit execution: onekgpd-individual-queries

Distinct increment over population-genetics-1000genomes: query individual genotypes/variants and select variant carriers or compute service-provided pairwise relatedness, beyond the retained IGSR population/sample discovery guide. Use scripts/onekgpd_meta.py for offline sample/population metadata from assets/kgpe.json. scripts/onekgpd_api.py needs the separately installed dnaerys client and queries the public db.dnaerys.org:443 service. Only public 1000 Genomes cohort identifiers are accepted; no patient data are needed. The MIT skill and metadata do not redistribute or relicense the proprietary Dnaerys database engine. Annotation sources/terms remain as described at https://dnaerys.org/online/ .

Install this complete folder with `dnakit skills install onekgpd-individual-queries --destination ./analysis-skills`. Inspect optional dependencies with `dnakit skills doctor onekgpd-individual-queries`. Resolve script/resource paths against the installed folder.

Import validation covers metadata, file integrity and Python syntax. It does not establish scientific accuracy, service availability or a successful pipeline run. Cite the upstream authors and the actual data/methods used; see UPSTREAM.md.
