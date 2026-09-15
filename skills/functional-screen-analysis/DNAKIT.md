# DNAKit integration

Canonical skill: `functional-screen-analysis`.

```bash
python -m pip install dnakit
dnakit skills doctor functional-screen-analysis
dnakit skills install functional-screen-analysis --agent codex
```

Follow SKILL.md for the scientific procedure and interpretation. DNAKit supplies the skill and its listed local resources. Real scientific libraries, command-line programs, models, reference databases and service credentials retain their own requirements.


## ENCODE

Use `dnakit.encode` for ENCODE queries, downloads and provenance tracking. The API retains the original options; `encode_search_experiments(...)` is now `dnakit.encode.search_experiments(...)`. It returns the original JSON text. Credential storage uses `dnakit[encode-credentials]`; public searches need no credentials.


## Provenance and validation

Original authors, license and revision remain in [UPSTREAM.md](UPSTREAM.md). DNAKit adaptation changes are recorded in the catalog. Packaging and interface checks do not establish scientific validity or completion of the full workflow.
