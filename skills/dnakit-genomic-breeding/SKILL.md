---
name: dnakit-genomic-breeding
description: Compute diploid genomic relationships and fixed-lambda RR-BLUP predictions with DNAKit. Use for genotype-to-phenotype breeding research when aligned dosage matrices and training phenotypes are available.
license: MIT
---

# Genomic relationships and breeding prediction

Install `pip install "dnakit[bio-analysis]"` when NumPy is missing.

1. Align samples, allele orientation and marker columns before computing. Input rows are individuals and columns are diploid allele dosages in [0, 2]; fractional dosages are allowed. Use null/NaN for missing values. The interface does not support polyploids.
2. Use `dnakit.genomics.genomic_relationship_matrix(genotypes)` for VanRaden method 1. The returned additive relationship matrix is not the kinship matrix (the conventional kinship is G/2).
3. For prediction, use `dnakit.genomics.genomic_prediction(train_genotypes, phenotypes, predict_genotypes, ridge_lambda=1.0)`. Choose lambda by training-only validation; 1.0 is an example, not a universal choice. Phenotypes must be finite and ordered with training rows.
4. Training allele means govern imputation and centering of prediction rows. All-missing and monomorphic training markers are excluded; inspect the returned metadata. Never tune parameters on the held-out targets.
5. Report prediction validation with a split appropriate to relatedness, population and breeding cycle. Fixed-lambda ridge regression does not estimate REML variance components, a multi-trait model, parental mating allocation, or a BOLE agent workflow.

Use the user's supplied organism, reference build and input files. Report inputs, parameters, backend versions, outputs and validation actually performed. Catalogue inclusion is not evidence that a dependency is installed or that an analysis has run.

Cite [VanRaden (2008)](https://doi.org/10.3168/jds.2007-0980), [Endelman (2011)](https://doi.org/10.3835/plantgenome2011.08.0024) and DNAKit's independent implementation. [BOLE, Wang et al. (2026)](https://doi.org/10.1016/j.jia.2026.07.011) is related workflow research; do not attribute DNAKit's implementation to unavailable BOLE source code.
