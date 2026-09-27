# DNAKit execution: kiro-healthomics-genomics

This is a transparent conversion of selected Kiro Markdown skill/steering guides, not an upstream SKILL.md or a bundled MCP implementation. Install the actual official server with `uvx awslabs.aws-healthomics-mcp-server`, configure an AWS profile/region, and inspect its live tools/list before mapping the upstream guide names. Reference: https://awslabs.github.io/mcp/servers/aws-healthomics-mcp-server/ . The AWS CLI exposes `aws omics create-workflow` and `aws omics start-run`. HealthOmics requires your IAM role, S3 input/output and container references; cloud resource creation/runs incur costs. No cloud writes were performed. Local life-sciences-* sample clients are not represented as validated scientific executors.

Install this complete folder with `dnakit skills install kiro-healthomics-genomics --destination ./analysis-skills`. Inspect optional dependencies with `dnakit skills doctor kiro-healthomics-genomics`. Resolve script/resource paths against the installed folder.

Import validation covers metadata, file integrity and Python syntax. It does not establish scientific accuracy, service availability or a successful pipeline run. Cite the upstream authors and the actual data/methods used; see UPSTREAM.md.
