---
name: dnakit-dna-nanotechnology
description: Construct an assigned scadnano DNA duplex and export caller-supplied scadnano geometry into oxDNA topology and configuration text using DNAKit.
license: MIT
---

# DNA nanostructure geometry

Install `pip install "dnakit[dna-nanotechnology]"`.

For a basic duplex use `design = dnakit.genomics.design_dna_duplex("AACCGGTT")`. Save the returned mapping as JSON for the scadnano editor when needed. For a supplied or generated design use `dnakit.genomics.export_oxdna_design(design)`; preserve both `configuration` and `topology` strings together.

Only linear strands are exported: the tested scadnano 0.21.1 exporter does not preserve circular closure bonds. All strands must have complete canonical sequences. Resolve unassigned or wildcard bases from the user's design intent before export; never assign arbitrary bases to make export pass. Inspect strand and nucleotide counts. The constructor only makes one straight duplex; it does not design a complete origami or optimize folding. Export prepares starting geometry and does not execute relaxation or oxDNA simulation. Molecular conclusions require an appropriate separately validated simulation and experimental context.

Use the user's supplied organism, reference build and input files. Report inputs, parameters, backend versions, outputs and validation actually performed. Catalogue inclusion is not evidence that a dependency is installed or that an analysis has run.

Cite [Doty, Lee and Stérin (2020), scadnano](https://doi.org/10.4230/LIPIcs.DNA.2020.9) and the [scadnano Python implementation](https://github.com/UC-Davis-molecular-computing/scadnano-python-package). If a subsequent oxDNA simulation is run, cite that engine and model as well.
