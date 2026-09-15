---
name: dnakit
description: Find and use curated DNA and genomics workflows from the DNAKit skill catalog. Use for DNA sequence processing, sequencing QC, alignment, assembly, variant analysis, GWAS, population genetics, epigenomics, CRISPR, and microbial genomics when the user needs an analysis workflow or wants to find an appropriate skill.
license: MIT
metadata:
  compatibility: Requires a DNAKit release containing the skill catalog, accessible through the dnakit CLI or DNAKit MCP server. Each selected workflow has additional dependencies.
---

## DNAKit execution

Use `dnakit skills doctor dnakit` to inspect scientific dependencies. Read [DNAKIT.md](DNAKIT.md) for this skill's DNAKit functions and bundled script entrypoints. Keep the workflow, input checks and reporting requirements below.


# DNAKit workflow discovery

## Find the workflow

1. Identify the task, input format, organism, reference build, and expected output from the user's request and available files. Ask only for information needed to proceed.
2. Search a few relevant task words with `dnakit skills search "variant annotation" --json --limit 5`. Chinese queries such as `变异注释` are supported. Use `dnakit skills categories --json` to narrow a broad request. With DNAKit MCP, use `search_dnakit_skills` instead.
3. Read the best match using `dnakit skills info SKILL_ID --json` and `dnakit skills read SKILL_ID`. With MCP, use `read_dnakit_skill`; follow `next_offset` when the response is truncated. Read the instructions before applying the workflow.

## Prepare and perform the task

- Check the selected skill's `runtime`, `required_skills`, provenance, and validation status. `dnakit skills doctor SKILL_ID` checks declared local dependencies only. It does not install dependencies or prove scientific correctness.
- To use supporting scripts or resources, install the complete selected folder with `dnakit skills install SKILL_ID --destination ./analysis-skills`. The command returns its exact path. Resolve that skill's resource paths from this directory. Text resources can also be read with `dnakit skills read SKILL_ID --file references/FILE.md` or the MCP reader's `resource` argument.
- Resolve cross-skill dependencies through the catalog's `required_skills` and search aliases. If a dependency is not bundled, use the recorded upstream source; do not invent a local file or silently substitute an unrelated workflow.
- Follow the selected workflow within the user's authorized task. Confirm reference build and coordinate conventions before joining genomic files. Check tool versions and actual command/API signatures when examples differ from the installed environment.
- Use tools available in the current environment. An upstream ToolUniverse, ClawBio, or other tool name is not automatically callable through DNAKit MCP. Install/configure the stated runtime when the task requires it, or report the specific missing dependency and complete independent work.
- Keep discovery selective: load only the chosen workflow and relevant references. Do not read the entire catalog into the conversation.

## Report the result

Record input files, reference build, parameters, software versions, output locations, and QC findings. Distinguish completed execution from a proposed or blocked workflow. Catalog inclusion is not execution validation.

Use `dnakit skills cite SKILL_ID` to obtain the pinned upstream software citation. Credit the original skill authors and the tools, models, datasets, and papers actually used. DNAKit provides curation, discovery, packaging, and integration; imported workflows retain their original licenses and authorship.
