# DNAKit integration

Canonical skill: `alphagenome-single-variant-analysis`.

```bash
python -m pip install dnakit
dnakit skills doctor alphagenome-single-variant-analysis
dnakit skills install alphagenome-single-variant-analysis --agent codex
```

Follow SKILL.md for the scientific procedure and interpretation. DNAKit supplies the skill and its listed local resources. Real scientific libraries, command-line programs, models, reference databases and service credentials retain their own requirements.


## Bundled scripts

Find scripts and run the relevant one with its documented arguments:

```bash
dnakit skills scripts alphagenome-single-variant-analysis
dnakit skills run alphagenome-single-variant-analysis --script scripts/analyze_ism.py -- <documented arguments>
```

The command resolves files from the installed DNAKit package; a checkout of the upstream skill collection is unnecessary. It forwards standard input/output and exit status. Python callers can use `dnakit.workflows.run(skill_id, script, arguments=[...], workdir=...)` for bounded captured output.

- `scripts/analyze_ism.py`
- `scripts/generate_ontology_mapping.py`
- `scripts/interpret_splicing.py`
- `scripts/lookup_gene_info.py`
- `scripts/resolve_ontology_terms.py`
- `scripts/visualize_genome_tracks.py`
- `scripts/visualize_variant_effects.py`

## Provenance and validation

Original authors, license and revision remain in [UPSTREAM.md](UPSTREAM.md). DNAKit adaptation changes are recorded in the catalog. Packaging and interface checks do not establish scientific validity or completion of the full workflow.
