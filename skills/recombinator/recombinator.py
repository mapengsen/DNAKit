#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
ClawBio Recombinator
Produce offspring genomes from parent pairs via meiotic recombination.

Usage:
    python recombinator.py --demo
    python recombinator.py --father einstein-g0 --mother anning-g0 --offspring 3
"""

import argparse
import json
import sys
from pathlib import Path

GENOMES_DIR = Path.cwd() / "genomebook-data" / "GENOMES"


def _import_recombinator():
    from dnakit.workflows.support import recombinator
    return recombinator


def load_genome(genome_id):
    """Load a genome by ID from the genomes directory."""
    path = GENOMES_DIR / f"{genome_id}.genome.json"
    if not path.exists():
        print(f"ERROR: Genome not found: {path}")
        sys.exit(1)
    with open(path) as f:
        return json.load(f)


def print_offspring(child):
    """Pretty-print an offspring genome summary."""
    print(f"\n{'=' * 60}")
    print(f"ID:     {child['id']}")
    print(f"Name:   {child['name']}")
    print(f"Sex:    {child['sex']} ({child['sex_chromosomes']})")
    print(f"Health: {child['health_score']}")
    print(f"Mutations: {len(child['mutations'])}")
    if child['mutations']:
        for m in child['mutations']:
            print(f"  - {m['locus']}: {m['from']}->{m['to']} ({m['type']}, from {m['parent']})")
    print(f"Conditions: {len(child['clinical_history'])}")
    for cond in child['clinical_history']:
        print(f"  - {cond['name']} ({cond['severity']}, fitness: {cond['fitness_cost']})")
    print(f"Top traits:")
    top = sorted(child['trait_scores'].items(), key=lambda x: x[1], reverse=True)[:5]
    for t, s in top:
        print(f"  - {t}: {s}")


def run_demo():
    """Run an explicitly synthetic one-locus example; no registry download is needed."""
    from dnakit.workflows import simulate_inheritance
    locus = {"chromosome": "1", "position": 1, "ref": "A", "alt": "C",
             "dominance": "additive", "effect_size": 1}
    father = {"id": "synthetic-father", "name": "Synthetic father", "sex": "Male",
              "loci": {"L1": {**locus, "alleles": ["A", "A"]}}}
    mother = {"id": "synthetic-mother", "name": "Synthetic mother", "sex": "Female",
              "loci": {"L1": {**locus, "alleles": ["C", "C"]}}}
    traits = {"traits": {"synthetic-allele-dosage": {"loci": [{"id": "L1", "ref": "A",
               "alt": "C", "dominance": "additive", "effect": 1}]}}}
    result = simulate_inheritance(father, mother, traits, {"diseases": {}},
                                  seed=42, mutation_rate=0, offspring_count=3)
    print("Synthetic per-locus inheritance demonstration; no real person or clinical prediction.")
    for child in result["offspring"]:
        print_offspring(child)
    print("[DEMO] Three synthetic offspring generated. No files written.")


def run_breed(father_id, mother_id, num_offspring, generation):
    """Breed a specific pair."""
    core = _import_recombinator()
    trait_reg, disease_reg = core.load_registries()

    father = load_genome(father_id)
    mother = load_genome(mother_id)

    print("=" * 70)
    print("RECOMBINATOR")
    print("=" * 70)
    print(f"Father: {father.get('name', father_id)} ({father_id})")
    print(f"Mother: {mother.get('name', mother_id)} ({mother_id})")
    print(f"Offspring: {num_offspring}, Generation: {generation}")

    children = core.breed_pair(
        father, mother, generation=generation, num_offspring=num_offspring,
        trait_reg=trait_reg, disease_reg=disease_reg
    )

    for child in children:
        print_offspring(child)

        # Write offspring genome
        out_path = GENOMES_DIR / f"{child['id']}.genome.json"
        with open(out_path, "w") as f:
            json.dump(child, f, indent=2)
        print(f"  -> Written: {out_path.name}")

    print(f"\n{len(children)} offspring written to {GENOMES_DIR}/")


def main():
    parser = argparse.ArgumentParser(
        description="Recombinator: Produce offspring via meiotic recombination"
    )
    parser.add_argument("--demo", action="store_true", help="Run a synthetic one-locus example")
    parser.add_argument("--father", type=str, help="Father genome ID (e.g. einstein-g0)")
    parser.add_argument("--mother", type=str, help="Mother genome ID (e.g. anning-g0)")
    parser.add_argument("--offspring", type=int, default=3, help="Number of offspring (default: 3)")
    parser.add_argument("--generation", type=int, default=1, help="Generation number (default: 1)")
    parser.add_argument("--data-dir", type=Path, default=Path.cwd() / "genomebook-data",
                        help="Directory containing GENOMES, trait_registry.json and disease_registry.json")
    parser.add_argument("--trait-registry", type=Path)
    parser.add_argument("--disease-registry", type=Path)
    args = parser.parse_args()
    global GENOMES_DIR
    GENOMES_DIR = args.data_dir / "GENOMES"
    core = _import_recombinator()
    core.TRAIT_REGISTRY = args.trait_registry or args.data_dir / "trait_registry.json"
    core.DISEASE_REGISTRY = args.disease_registry or args.data_dir / "disease_registry.json"

    if args.demo:
        run_demo()
    elif args.father and args.mother:
        run_breed(args.father, args.mother, args.offspring, args.generation)
    else:
        parser.print_help()
        print("\nUse --demo for a quick demo, or specify --father and --mother.")


if __name__ == "__main__":
    main()
