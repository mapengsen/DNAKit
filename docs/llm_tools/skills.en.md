# DNAKit Skills: search, install, and use

DNAKit collects DNA/genomics workflows with offline discovery, selective installation, provenance, and MCP access. A Skill contains instructions and supporting resources. Its analysis software, reference data, and APIs remain separate dependencies.

Start with the **`dnakit` discovery skill**. Describe your task; the Agent searches the catalog, reads a suitable workflow, and uses the available tools to perform it.

## Install DNAKit and the discovery skill

This feature requires a DNAKit release containing the skill catalog. From the current source checkout:

```bash
python -m pip install -e "."
dnakit skills stats
```

After a release containing this feature is published, users can obtain it with `python -m pip install -U dnakit`. Discovery and installation do not require a model API or MCP.

In your analysis project, choose one client:

```bash
dnakit skills install dnakit --agent codex
# Or:
dnakit skills install dnakit --agent claude-code
# Or:
dnakit skills install dnakit --agent cursor
```

| Client | Project directory | User directory with `--global` |
| --- | --- | --- |
| Codex | `.agents/skills/dnakit/` | `~/.agents/skills/dnakit/` |
| Claude Code | `.claude/skills/dnakit/` | `~/.claude/skills/dnakit/` |
| Cursor | `.cursor/skills/dnakit/` | `~/.cursor/skills/dnakit/` |

See the official [Codex](https://learn.chatgpt.com/docs/build-skills), [Claude Code](https://code.claude.com/docs/en/skills), and [Cursor](https://cursor.com/docs/skills) documentation. Start a new client session after installation, and make the DNAKit Python environment accessible to the client.

For example, ask Codex:

```text
$dnakit Convert BED interval chr1 999 1000 to GFF coordinates.
Find a suitable skill, check dependencies, run the conversion, and explain the conventions.
```

The expected route is discovery of `genomic-coordinates`, reading its instructions, preparing its runtime, and running its conversion script. The result is `chr1 1000 1000`, a 1 bp interval.

## Find and install an individual workflow

```bash
dnakit skills search "variant annotation" --limit 5
dnakit skills search "变异注释" --limit 5
dnakit skills search "GWAS" --category gwas --limit 5
dnakit skills categories
dnakit skills info variant-annotation
dnakit skills doctor variant-annotation
dnakit skills install variant-annotation --agent codex
```

Search uses local keyword matching and curated Chinese query aliases. It makes no model API calls; summaries and workflow text returned to an Agent still consume conversation tokens. Skill discovery finds workflow instructions; tool discovery finds callable functions.

`doctor` probes declared Python modules, commands, and environment-variable presence without printing secret values. Read the workflow for complete version, R package, reference-data, API, and model/data-license requirements. Passing local checks is not execution validation.

Installation includes the selected skill's supporting files and licenses. Existing unmanaged directories are protected. Identical installations are unchanged; replacing a modified DNAKit installation requires explicit `--force`. Analysis dependencies are not installed by this command.

## Run an offline example

The upstream `genomic-coordinates` skill declares Python 3.11+:

```bash
dnakit skills install genomic-coordinates --destination ./analysis-skills
python3 analysis-skills/genomic-coordinates/scripts/convert_coords.py \
  --from bed --to gff chr1 999 1000 --format json
```

The output has `start=1000`, `end=1000`, and `length=1`: BED is 0-based half-open, while GFF is 1-based inclusive. Converting conventions is distinct from mapping between reference genome builds.

## Use MCP

After following the [MCP tutorial](mcp.md), an Agent can call:

```text
search_dnakit_skills(query="variant annotation", limit=5)
read_dnakit_skill(skill_id="variant-annotation")
```

Use a relative path from the response's `files` list as `resource` to read a supporting document. Continue with `next_offset` when provided. Reading text through MCP does not install scripts on the client; install the selected folder in the execution environment when scripts are needed.

The full catalog now has DNAKit integration guides and maintained adaptations.
Collection adapters and bundled scripts use the [unified DNAKit APIs](native_skills.md),
while actual scientific software and models retain their dependencies.
The `sequence-retrieval` workflow uses native NCBI/ENA adapters. FASTA, NCBI GenBank and ENA EMBL use the base package;
cross-format conversion requires `dnakit[bio-formats]`. See
[sequence search](../api/features/22_database_query.md) and
[sequence downloads](../api/features/15_download.md#accession-sequences).

```bash
dnakit skills search "sequence retrieval" --limit 5
dnakit skills install sequence-retrieval --agent codex
dnakit skills read sequence-retrieval --file examples.md
```

The old `tooluniverse-sequence-retrieval` and `public-db-lookup` identifiers remain
aliases. Existing client installation folders are not automatically moved or
deleted; remove obsolete copies when appropriate to avoid duplicate discovery.
Use each skill's documentation and `skills info` runtime metadata to check its status.

## Provenance, deduplication, and citation

This snapshot imports workflows from bioSkills, ToolUniverse, ClawBio, SciAgent, K-Dense, and NVIDIA BioNeMo, alongside DNAKit's discovery skill. Reviewed task/stage equivalence rules choose one implementation per duplicate group. Different input conditions, composite workflows, and their component stages can remain distinct. Registered alternative names resolve to their retained implementation.

```bash
dnakit skills stats
dnakit skills list --source tooluniverse --limit 20
dnakit skills info genomic-coordinates --json
dnakit skills cite genomic-coordinates
dnakit skills cite genomic-coordinates --format bibtex
```

Each upstream skill includes `UPSTREAM.md`, license notices, a pinned Git revision, and resource checksums. See [`registry/`](https://github.com/mapengsen/DNAKit/tree/main/registry) for selection decisions and [`registry/skills.csv`](https://github.com/mapengsen/DNAKit/blob/main/registry/skills.csv) for a publication-friendly inventory.

DNAKit code remains MIT; bundled third-party content retains its MIT, Apache-2.0, CC-BY-4.0, or AGPL-3.0-only terms. Independent ENCODE skills remain AGPL with complete editable sources. DeepMind software and explanatory content retain Apache-2.0 and CC-BY-4.0 respectively. External MCP services, model weights, and databases require their own setup and terms.

The expanded collection includes DeepMind science-skills, PKU-EMBL Metagenomics-Skills, omics-skills, ENCODE Toolkit, and other original sources. It covers DNA models, microbial and viral genomes, plasmids, mobile genetic elements, and regulatory analysis. Search, for example, with `dnakit skills search "DNABERT" --category models` or `dnakit skills search "integron"`.

See [`registry/discovery.md`](https://github.com/mapengsen/DNAKit/blob/main/registry/discovery.md) for search terms, pinned sources, and exclusions. Sources without an established content license, current noncommercial or competing-use restrictions, or missing essential scripts were not copied into this general-use distribution. This is a reproducible public-web snapshot, not a claim of exhaustive internet coverage.

Credit DNAKit's curation, discovery, packaging, and integration separately from original workflow authorship. Cite upstream software, methods, models, and datasets actually used. Inclusion and checksum validation do not establish scientific correctness. The catalog and `registry/validation.json` record the scope of execution checks.

## Unified execution

Every bundled skill now includes `DNAKIT.md` with its DNAKit APIs, packaged script entrypoints and retained scientific dependencies. Use `dnakit skills scripts ID` and `dnakit skills run ID --script FILE -- <arguments>` without an upstream checkout. See the [native interface guide](native_skills.md).
