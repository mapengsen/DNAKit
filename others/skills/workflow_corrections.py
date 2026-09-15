"""Reviewed packaging corrections for composite workflow scripts and demos."""

from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
TEMPLATES = ROOT / "others/skills/adaptations"


def main():
    path = TEMPLATES / "nfcore-sarek-wrapper/nfcore_sarek_wrapper.py"
    text = path.read_text()
    for name in (
        "clinical-variant-reporter",
        "clinical-trial-finder",
        "omics-target-evidence-mapper",
    ):
        original = f"python3 skills/{name}/{name.replace('-', '_')}.py "
        text = text.replace(original, f"dnakit skills support-run {name} -- ")
    # The optional downstream stages consume different inputs; retain their actual
    # CLI contracts instead of forwarding a directory of annotated VCFs to all of them.
    text = text.replace(
        '"dnakit skills support-run clinical-trial-finder -- "\n                f"--input {shlex.quote(str(annotation_dir))} --output <dir>"',
        '"dnakit skills support-run clinical-trial-finder -- "\n                "--gene <gene_symbol> --output <dir>"',
    )
    text = text.replace(
        '"dnakit skills support-run omics-target-evidence-mapper -- "\n                f"--input {shlex.quote(str(annotation_dir))} --output <dir>"',
        '"dnakit skills support-run omics-target-evidence-mapper -- "\n                "--gene <gene_symbol> --output <dir>"',
    )
    for language in ("en", "es"):
        old = (
            f'"python3 skills/wes-clinical-report-{language}/wes_clinical_report_{language}.py "\n'
            '                f"--report-dir {shlex.quote(str(output_dir))} --output-dir <pdf_dir>"'
        )
        new = f'"dnakit skills report --language {language} "\n                "--input <wes_report.md> --output <pdf_dir>"'
        text = text.replace(old, new)
    path.write_text(text)
    path = TEMPLATES / "nfcore-sarek-wrapper/reporting.py"
    text = path.read_text().replace(
        "# Requires a ClawBio checkout. Uncomment, set CLAWBIO_REPO and an output dir:",
        "# Requires DNAKit. Uncomment and select an output directory:",
    )
    text = text.replace(
        'CLAWBIO_REPO=/path/to/ClawBio python3 "$CLAWBIO_REPO/clawbio.py" run sarek-pipeline',
        "dnakit skills run nfcore-sarek-wrapper --script nfcore_sarek_wrapper.py --",
    )
    path.write_text(text)
    path = TEMPLATES / "illumina-bridge/illumina_bridge.py"
    text = path.read_text()
    for old, new in {
        "python clawbio.py run equity": "dnakit skills run equity-scorer --script equity_scorer.py --",
        "python clawbio.py run clinpgx": "dnakit skills run clinpgx --script clinpgx.py --",
        "python clawbio.py run gwas --rsid <rsid> --output <dir>": 'dnakit database call gwas_get_associations_for_snp \'{\\"rs_id\\": \\"<rsid>\\"}\'',
        "python clawbio.py run {SKILL_ALIAS}": "dnakit skills run illumina-bridge --script illumina_bridge.py --",
    }.items():
        text = text.replace(old, new)
    path.write_text(text)
    path = TEMPLATES / "recombinator/recombinator.py"
    text = (ROOT / "skills/recombinator/recombinator.py").read_text()
    if "_PROJECT_ROOT =" not in text:
        return  # The public skill already uses the maintained DNAKit engine template.
    a, b = text.index("_PROJECT_ROOT ="), text.index("def load_genome(")
    text = (
        text[:a]
        + """GENOMES_DIR = Path.cwd() / "genomebook-data" / "GENOMES"


def _import_recombinator():
    from dnakit.workflows.support import recombinator
    return recombinator


"""
        + text[b:]
    )
    a, b = text.index("def run_demo():"), text.index("def run_breed(")
    text = (
        text[:a]
        + '''def run_demo():
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


'''
        + text[b:]
    )
    text = text.replace(
        'help="Demo with Einstein x Anning"', 'help="Run a synthetic one-locus example"'
    )
    text = text.replace(
        "    args = parser.parse_args()",
        '''    parser.add_argument("--data-dir", type=Path, default=Path.cwd() / "genomebook-data",
                        help="Directory containing GENOMES, trait_registry.json and disease_registry.json")
    parser.add_argument("--trait-registry", type=Path)
    parser.add_argument("--disease-registry", type=Path)
    args = parser.parse_args()
    global GENOMES_DIR
    GENOMES_DIR = args.data_dir / "GENOMES"
    core = _import_recombinator()
    core.TRAIT_REGISTRY = args.trait_registry or args.data_dir / "trait_registry.json"
    core.DISEASE_REGISTRY = args.disease_registry or args.data_dir / "disease_registry.json"''',
    )
    path.write_text(text)
    guide = TEMPLATES / "recombinator/DNAKIT.md"
    text = guide.read_text()
    if "## Simulation data" not in text:
        guide.write_text(
            text + "\n## Simulation data\n\nDNAKit bundles the attributed inheritance engine. "
            "Use `--demo` for a synthetic one-locus example. For supplied parent genomes, use `--data-dir` "
            "with `GENOMES/`, `trait_registry.json` and `disease_registry.json`, or explicit registry paths. "
            "The original historical-person demo genomes were absent from the pinned repository and are not fabricated here. "
            "The Python API is `dnakit.workflows.simulate_inheritance(...)`, with optional offspring count, mutation rate and seed. "
            "This is independent allele inheritance, not chromosome-scale crossover or a clinical predictor.\n"
        )


if __name__ == "__main__":
    main()
