# DNAKit execution: okg-search

This skill is curated from `svdorn/statgen-skills` at `dc5f0c60d2995aebfe58b3f8a07b46f27562f4e0`. DNAKit provides discovery,
packaging, dependency inspection and provenance; the scientific backend remains upstream.
Read [UPSTREAM.md](UPSTREAM.md) and [the original project README](references/upstream-readme.md).

Inspect dependencies with `dnakit skills doctor okg-search` after catalog integration.
Install the complete folder with `dnakit skills install okg-search --destination ./analysis-skills`.
Resolve scripts/resources against that installed folder; run `python scripts/NAME.py --help`
from the folder before using a script. Do not assume an upstream MCP tool is exposed by DNAKit.
The import has metadata, resource-hash and Python-syntax checks. Scientific execution and
prediction accuracy have NOT been validated by importing this workflow.

Run the bundled script from this skill's directory in a dedicated analysis
environment. Upstream helpers may install Python/R packages, clone scientific tools or
download large reference panels on first ANALYSIS use; importing this skill performs none
of these operations. Inspect arguments and reference paths before an analysis run.

OKG integration is a separate service: when using it, supply your own `OKG_REPO`
(containing deployments/statgen-analysis/server.py) and `OKG_DSN`. The repository's
example localhost DSN is not a supplied database. The OKG schema/data/server are not
bundled here. `okg-search` and the current PRS/TWAS metadata paths depend on this
configuration; other workflows have standalone paths documented upstream.
R/PLINK/LDSC/SBayesRC/FUSION/SuSiE and LD/GTEx/reference files are separate dependencies,
not implementations supplied by the skill package. Preserve build, allele and ancestry
matching; report actual solver/reference versions. Runtime downloads can be unpinned
upstream even though this skill is pinned, so record resolved versions for reproducibility.
