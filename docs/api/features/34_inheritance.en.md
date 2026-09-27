# Inheritance simulation

<!-- BEGIN GENERATED API: inheritance -->

<span id="inheritance"></span>

## 1) Simulate allele inheritance at independent loci {#dnakit.workflows.simulate_inheritance}

- **Purpose:** Simulate allele inheritance at independent loci.
- **Method:** Sample one allele from each parent at each independent locus, then apply the specified mutation probability.
- **API：** `dnakit.workflows.simulate_inheritance(father, mother, trait_registry, disease_registry, *, offspring_count=3, generation=1, mutation_rate=0.001, seed=None, verbose=True)`
- **Input:** Provide `father`, `mother`, `trait_registry`, `disease_registry`. This example uses an explicit synthetic locus with mutation disabled and no remote data.
- **Output:** Mapping with offspring, seed, mutation_rate, scope and upstream; independent-locus simulation, not chromosome crossover or clinical prediction.
- **Example code:**

```python
import dnakit

example_father = {'id': 'father',
 'sex': 'Male',
 'loci': {'L1': {'chromosome': '1',
                 'position': 1,
                 'ref': 'A',
                 'alt': 'C',
                 'dominance': 'additive',
                 'effect_size': 1,
                 'alleles': ['A', 'A']}}}

example_mother = {'id': 'mother',
 'sex': 'Female',
 'loci': {'L1': {'chromosome': '1',
                 'position': 1,
                 'ref': 'A',
                 'alt': 'C',
                 'dominance': 'additive',
                 'effect_size': 1,
                 'alleles': ['C', 'C']}}}

example_trait_registry = {'traits': {'dosage': {'loci': [{'chromosome': '1',
                                 'position': 1,
                                 'ref': 'A',
                                 'alt': 'C',
                                 'dominance': 'additive',
                                 'effect_size': 1,
                                 'id': 'L1',
                                 'effect': 1}]}}}

result = dnakit.workflows.simulate_inheritance(
    father=example_father,
    mother=example_mother,
    trait_registry=example_trait_registry,
    disease_registry={'diseases': {}},
    offspring_count=2,
    mutation_rate=0,
    seed=13,
    verbose=False,
)
print(result)
```

- **Example output:** Mapping with offspring, seed, mutation_rate, scope and upstream; independent-locus simulation, not chromosome crossover or clinical prediction.

Response-shape excerpt; `<…>` denotes type placeholders, not a captured query result.

```json
{
  "offspring": "<list>",
  "seed": "<int | None>",
  "mutation_rate": "<float>",
  "scope": "<str>",
  "upstream": "<str>"
}
```

[Full parameters and usage notes](../../faq/api/34_inheritance.md#dnakit.workflows.simulate_inheritance-parameters)

<!-- END GENERATED API: inheritance -->
