# DNAKit execution: tangermeme

This skill is curated from `jmschrei/tangermeme` at `d5c51b83d374576b64effb9284fc79a71a604734`. DNAKit provides discovery,
packaging, dependency inspection and provenance; the scientific backend remains upstream.
Read [UPSTREAM.md](UPSTREAM.md) and [the original project README](references/upstream-readme.md).

Inspect dependencies with `dnakit skills doctor tangermeme` after catalog integration.
Install the complete folder with `dnakit skills install tangermeme --destination ./analysis-skills`.
Resolve scripts/resources against that installed folder; run `python scripts/NAME.py --help`
from the folder before using a script. Do not assume an upstream MCP tool is exposed by DNAKit.
The import has metadata, resource-hash and Python-syntax checks. Scientific execution and
prediction accuracy have NOT been validated by importing this workflow.

Install `tangermeme`/PyTorch in a separate model environment.
The source snapshot
can be installed with `python -m pip install git+https://github.com/jmschrei/tangermeme.git@d5c51b83d374576b64effb9284fc79a71a604734`.
Use the real `tangermeme.predict`, `tangermeme.deep_lift_shap`,
`tangermeme.variant_effect`, and `tangermeme.design` modules as documented in references.
Supply a compatible trained model and its preprocessing; no pretrained model or weights
are included here. GPU, memelite and model-specific dependencies are optional by task.
