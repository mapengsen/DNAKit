# Inheritance simulation: parameters and notes

See the [tutorial](../../api/features/34_inheritance.md) for usage and examples.

Each function includes its purpose, API, inputs, outputs, example code and result description. Prepare any required credentials, local files and optional dependencies; remote results vary with database releases and query time.

## Simulate allele inheritance at independent loci {#dnakit.workflows.simulate_inheritance-parameters}

`dnakit.workflows.simulate_inheritance(father, mother, trait_registry, disease_registry, *, offspring_count=3, generation=1, mutation_rate=0.001, seed=None, verbose=True)`

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `father` | dict&#91;str, Any&#93; | Required | Father record with sex=Male and Genomebook genotype fields. |
| `mother` | dict&#91;str, Any&#93; | Required | Mother record with sex=Female and Genomebook genotype fields. |
| `trait_registry` | dict&#91;str, Any&#93; | Required | Explicit trait and locus registry. |
| `disease_registry` | dict&#91;str, Any&#93; | Required | Explicit disease locus registry. |
| `offspring_count` | int | 3 | Number of simulated offspring, 1–10000. |
| `generation` | int | 1 | Positive integer generation label. |
| `mutation_rate` | float | 0.001 | Per-locus simulated mutation rate, between 0 and 1. |
| `seed` | int &#124; None | None | Local random seed; does not modify global RNG state. |
| `verbose` | bool | True | Show DNAKit progress messages without changing results. |
