# DNAKit execution: qtl-candidate-gene

This skill is curated from `tayloranthonyanderson/haploqtl` at `210eace2912c7bf11d119cb503c89f8f6cad479f`. DNAKit provides discovery,
packaging, dependency inspection and provenance; the scientific backend remains upstream.
Read [UPSTREAM.md](UPSTREAM.md) and [the original project README](references/upstream-readme.md).

Inspect dependencies with `dnakit skills doctor qtl-candidate-gene` after catalog integration.
Install the complete folder with `dnakit skills install qtl-candidate-gene --destination ./analysis-skills`.
Resolve scripts/resources against that installed folder; run `python scripts/NAME.py --help`
from the folder before using a script. Do not assume an upstream MCP tool is exposed by DNAKit.
The import has metadata, resource-hash and Python-syntax checks. Scientific execution and
prediction accuracy have NOT been validated by importing this workflow.

The full upstream script/reference/data folder is retained. Interval lookup uses
the bundled ITAG4.1 slice offline. Supply your own matching full GFF outside EB-5/EB-9.
MAS analysis additionally requires NumPy/scikit-allel and an input VCF with stated groups;
the disease-group defaults are only the published EB-9 example. UniProt/PubMed steps
use network; optional citation-support checking needs anthropic and ANTHROPIC_API_KEY.
The core interval scripts do not require installing the whole haploqtl project.
Cite Anderson et al. (2024), DOI 10.1111/tpj.16495, and the SGN/UniProt/NCBI data used.
See the retained [CITATION.cff](CITATION.cff).
