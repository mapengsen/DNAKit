# DNAKit execution: gigwa-explore-instance

This skill is curated from `gkanogiannis/Gigwa-MCP` at `d6ccdd1a5eb19b8a0cc224b96ad459868d9c4aaf`. DNAKit provides discovery,
packaging, dependency inspection and provenance; the scientific backend remains upstream.
Read [UPSTREAM.md](UPSTREAM.md) and [the original project README](references/upstream-readme.md).

Inspect dependencies with `dnakit skills doctor gigwa-explore-instance` after catalog integration.
Install the complete folder with `dnakit skills install gigwa-explore-instance --destination ./analysis-skills`.
Resolve scripts/resources against that installed folder; run `python scripts/NAME.py --help`
from the folder before using a script. Do not assume an upstream MCP tool is exposed by DNAKit.
The import has metadata, resource-hash and Python-syntax checks. Scientific execution and
prediction accuracy have NOT been validated by importing this workflow.

Install `gigwa-mcp` separately, connect its MCP server, and configure a local or
remote Gigwa instance. `GIGWA_URL` selects the instance; credentials are optional for
public reads and required for private data/imports. Inputs must be readable by the server.
Import only into the dataset/server selected by the user. `DNAKit` does not install Gigwa.
Read [NOTICE](NOTICE) for the distinction between Gigwa and this MCP wrapper.

Scientific interpretation: VanRaden G is the additive relationship matrix, whereas
coancestry is conventionally G/2. No universal G threshold proves sibling/clone status.
Complete call rate alone does not establish a bad import. UPGMA from IBS is a sample
distance dendrogram; PCA plus K-means is a structure heuristic, not ADMIXTURE inference.
Apply these qualifications when reading the upstream example wording.
