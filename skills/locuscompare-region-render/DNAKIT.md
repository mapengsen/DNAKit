# DNAKit integration

Canonical skill: `locuscompare-region-render`.

```bash
python -m pip install dnakit
dnakit skills doctor locuscompare-region-render
dnakit skills install locuscompare-region-render --agent codex
```

Follow SKILL.md for the scientific procedure and interpretation. DNAKit supplies the skill and its listed local resources. Real scientific libraries, command-line programs, models, reference databases and service credentials retain their own requirements.


## Bundled scripts

Find scripts and run the relevant one with its documented arguments:

```bash
dnakit skills scripts locuscompare-region-render
dnakit skills run locuscompare-region-render --script _fetchers/__init__.py -- <documented arguments>
```

The command resolves files from the installed DNAKit package; a checkout of the upstream skill collection is unnecessary. It forwards standard input/output and exit status. Python callers can use `dnakit.workflows.run(skill_id, script, arguments=[...], workdir=...)` for bounded captured output.

- `_fetchers/__init__.py`
- `_fetchers/gencode_ondemand.py`
- `_prefetched.py`
- `_wald_ratio_types.py`
- `cli.py`
- `commands.sh`
- `examples/01_synthetic_demo/generate_synthetic_fixtures.py`
- `examples/recipes/finngen_direct/harmonise.sh`
- `examples/recipes/gtex_v10_direct/harmonise.sh`
- `examples/recipes/pan_ukbb_direct/harmonise.sh`
- `examples/recipes/ukb_ppp_pqtl/harmonise.sh`
- `locuscompare_region_render.py`
- `regional_plot.py`

## Provenance and validation

Original authors, license and revision remain in [UPSTREAM.md](UPSTREAM.md). DNAKit adaptation changes are recorded in the catalog. Packaging and interface checks do not establish scientific validity or completion of the full workflow.
