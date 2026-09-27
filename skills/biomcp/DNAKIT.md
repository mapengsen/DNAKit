# DNAKit execution: biomcp

This skill is curated from `genomoncology/biomcp` at `3b1d6452ab63741562d2ef97495a97281b304729`. DNAKit provides discovery,
packaging, dependency inspection and provenance; the scientific backend remains upstream.
Read [UPSTREAM.md](UPSTREAM.md) and [the original project README](references/upstream-readme.md).

Inspect dependencies with `dnakit skills doctor biomcp` after catalog integration.
Install the complete folder with `dnakit skills install biomcp --destination ./analysis-skills`.
Resolve scripts/resources against that installed folder; run `python scripts/NAME.py --help`
from the folder before using a script. Do not assume an upstream MCP tool is exposed by DNAKit.
The import has metadata, resource-hash and Python-syntax checks. Scientific execution and
prediction accuracy have NOT been validated by importing this workflow.

Install the upstream BioMCP CLI separately according to its README, then inspect
`biomcp list` and `biomcp skill list`. The pinned project may use a different runtime
from older Python packages with a similar name; do not substitute an unverified package.
Gene/variant/PGx/literature functions query upstream services and may use downloaded caches
or credentials for particular providers. Use the genomic/biomedical research playbooks
actually relevant to the request; source data, evidence and model limitations remain visible.
SDLC issue-filing instructions are excluded. Cite the retained CITATION.cff and data sources.
