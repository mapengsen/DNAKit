---
name: dnakit-prokaryotic-genes
description: Predict prokaryotic and metagenomic CDS with optional Pyrodigal through DNAKit, returning coordinates, coding sequences and proteins. Use for bacterial or archaeal sequence annotation.
license: MIT
---

# Prokaryotic gene prediction

Install `pip install "dnakit[genomics-genes]"`.

Use `dnakit.genomics.predict_prokaryotic_genes(sequence, sequence_id="contig", meta=True)` and call `.to_dict()` for JSON. Metagenomic mode uses pretrained models. For single-genome training use `meta=False` and a suitable training sequence (roughly 100 kb or longer); inspect Pyrodigal's training requirements before selecting that mode. Explicit translation tables apply only where supported by that mode.

- Supply one linear contig per call and retain its identifier. Returned coordinates are 0-based half-open on the input sequence; reverse-strand CDS is reverse-complemented. Do not reinterpret these as 1-based GFF coordinates.
- Ambiguous IUPAC bases become N; report the normalization metadata and partial genes.
- This is prokaryotic CDS prediction. Eukaryotic exons, introns, splice forms and functional gene names require distinct workflows. Circular wrapping genes require a deliberately chosen linearization.

Use the user's supplied organism, reference build and input files. Report inputs, parameters, backend versions, outputs and validation actually performed. Catalogue inclusion is not evidence that a dependency is installed or that an analysis has run.

Cite Martin Larralde's [Pyrodigal](https://doi.org/10.21105/joss.04296) and Hyatt et al.'s [Prodigal](https://doi.org/10.1186/1471-2105-11-119). Pyrodigal is an optional separately licensed dependency, not vendored source.
