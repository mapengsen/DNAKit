---
name: dnakit-methylation-summary
description: Read modkit extended bedMethyl and compute coverage-weighted methylation summaries with DNAKit. Use for nanopore modification counts, including separate 5mC and 5hmC groups.
license: MIT
---

# Modification counts and methylation summaries

Use `dnakit.genomics.read_bedmethyl(path)` for modkit's 18-column extended bedMethyl; gzip/BGZF input is accepted. The reader bounds record count, decompressed bytes and line length. If a real dataset exceeds a bound, select a deliberate new bound or partition the input.

Pass the records to `dnakit.genomics.summarize_methylation(records, min_coverage=10)`; choose the coverage threshold for the experiment. Preserve chromosome, strand, modification code and motif strata. The reported fraction is sum(modified counts)/sum(valid coverage), never an unweighted mean of rounded per-site percentages. Zero coverage yields no fraction.

Column 11 is a percentage (0–100), despite some abbreviated upstream prose calling it a fraction. Modified, canonical and other modified counts must sum to valid coverage. Use the interface's validation errors to identify malformed input; do not silently discard inconsistent records. Combined-strand rows cannot be split back into independent strand measurements.

Generating bedMethyl from modBAM requires the separate `modkit-pileup` external analysis. Descriptive summaries do not perform DMR testing, infer causality, or predict disease.

Use the user's supplied organism, reference build and input files. Report inputs, parameters, backend versions, outputs and validation actually performed. Catalogue inclusion is not evidence that a dependency is installed or that an analysis has run.

Cite the [Oxford Nanopore modkit authors and format documentation](https://nanoporetech.github.io/modkit/intro_pileup.html) and DNAKit's independently written parser. Keep modkit's own software terms separate from this skill's MIT license.
