# Upstream provenance: sequence-retrieval

- Source: [mims-harvard/ToolUniverse](https://github.com/mims-harvard/ToolUniverse/tree/7a0ceb2f3a2ae28594e55ec27344e839df8bc9f7/skills/tooluniverse-sequence-retrieval)
- Revision: `7a0ceb2f3a2ae28594e55ec27344e839df8bc9f7`
- Original name: `tooluniverse-sequence-retrieval`
- Original path: `skills/tooluniverse-sequence-retrieval`
- Attribution: ToolUniverse contributors, mims-harvard
- Content license: **Apache-2.0**. Repository content license applies; underlying tools/data are not bundled.
- Original import snapshot: 2026-09-12
- DNAKit adaptation: 2026-09-13. Modified `SKILL.md`, `examples.md` and `CHECKLIST.md` to use DNAKit APIs, clarify accession routing and biological evidence, and implement native/download-conversion format choices. The full retrieval, disambiguation and reporting workflow is retained.
- Identifier renamed to `sequence-retrieval`; the original identifier and `public-db-lookup` remain catalog aliases.
- DNAKit additions: this provenance notice, copies of source notices in `licenses/`, and native Python adapters in the DNAKit package.
- Original metadata, upstream file hashes and source revision remain in `registry/catalog.json`. Reviewed adaptation templates are maintained under `others/skills/adaptations/sequence-retrieval/`.

## Runtime and validation

The adapted workflow uses DNAKit's NCBI and ENA HTTP adapters; it does not import
ToolUniverse. Public services require network access. Biopython remains an optional
dependency for annotated-format conversion only (`dnakit[bio-formats]`). Database,
library and data terms remain separate; this bundle contains no database mirror.

The original workflow was not executed during import. Current execution evidence
must be assessed separately from metadata/checksum validation; inspect
`dnakit skills info sequence-retrieval --json` and DNAKit's adapter tests. Passing
adapter checks does not establish Agent task success or scientific conclusions.

## Citation

Cite mims-harvard/ToolUniverse at the pinned revision and the underlying methods
actually used. Consult the upstream README and CITATION.cff, when provided, for
the project's preferred paper citation.

ToolUniverse contributors, mims-harvard. **tooluniverse-sequence-retrieval**.
Revision `7a0ceb2f3a2ae28594e55ec27344e839df8bc9f7`.
[Original Skill](https://github.com/mims-harvard/ToolUniverse/tree/7a0ceb2f3a2ae28594e55ec27344e839df8bc9f7/skills/tooluniverse-sequence-retrieval).
Adapted as **sequence-retrieval** by DNAKit, 2026-09-13.

Also cite the DNAKit version, NCBI/ENA records and versions, and any libraries or
papers used in the actual analysis. Inclusion does not imply upstream endorsement.


## DNAKit adaptation — 2026-09-13

Distributed as `sequence-retrieval`; previous identifier `sequence-retrieval` remains an alias when different. DNAKit changed discovery links, package calls and execution guidance. Original authorship, licenses, citations and original file hashes are retained. Scientific software/models remain their original implementations. See DNAKIT.md and registry/adaptations.json.
