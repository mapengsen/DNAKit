# Microbiome

This category contains **21** documented candidate capabilities; 12 compute real results natively (abundance-table handling, diversity, k-mer ANI estimate, co-occurrence network, differential abundance, dereplication, plotting and export) and 9 always return `status: "not_implemented"` (`assemble`, `bin`, `build_microbial_gene_catalog`, `build_pathogen_transmission_tree`, `calculate_digital_ddh`, `classify`, `infer_viral_quasispecies`, `predict`, `track`) because they need an external program, a reference database or trained models. `calculate_differential_abundance` additionally returns `not_implemented` when scipy is not installed. Two results are explicit heuristics: `calculate_ani` (and `dereplicate` in `ani` mode) is a Mash-formula estimate from exact k-mer sets, not an alignment-based ANI, and the co-occurrence network uses plain correlations of relative abundance. The documented functions are exported through the formal DNAKit namespace; functions with a compatible JSON call plan are also exposed through the Agent/MCP registry.

## 1) annotate

- **Purpose:** Parse taxonomy lineage strings into ranked fields (domain to species) and report the lowest assigned rank.

- **Method:** Pure text parsing: a lineage is split at `;` or `|`. If every level carries a rank prefix (`k__`/`d__`, `p__`, `c__`, `o__`, `f__`, `g__`, `s__`; `k__` and `d__` both map to `domain`) the prefix decides the rank; otherwise the levels are read positionally as domain, phylum, class, order, family, genus, species (more than 7 unprefixed levels is an error). Empty or placeholder names (`__`, `Unassigned`, `unclassified`, `unknown`, `uncultured`, `ambiguous`, `incertae sedis`, ...) become `null`. This does not classify sequences; it only interprets lineage text that is already assigned.

- **API:** `dnakit.microbiome.annotate(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input`: one lineage string, a list of lineage strings (ids `lineage1`, ...) or a mapping `{feature_id: lineage}`, at most 1,000,000 lineages. `options`: `entity` (only `"taxonomy"`, the default).

- **Output:** The envelope carries `status` (`"ok"`), `module` (`"microbiome"`), `operation`, `provenance` (`implementation`, `backend`), `input_summary` (a generic text-level summary: record count, text length, base counts, SHA-256; a mapping input counts as one record), `warnings` and `result`. `result` is a list with one entry per lineage: `id`, `lineage`, `ranks` (a mapping from `domain`, `phylum`, `class`, `order`, `family`, `genus`, `species` to a name or `null`), `lowest_rank` and `lowest_taxon`.

- **Example code:**

```python
import json
from dnakit.microbiome import annotate

lineages = {
    "otu1": "k__Bacteria;p__Firmicutes;c__Bacilli;o__Lactobacillales;f__Lactobacillaceae;g__Lactobacillus;s__",
    "otu2": "Archaea;Euryarchaeota;Methanobacteria",
}
result = annotate(lineages)
print(json.dumps(result, indent=2))
```

- **Example output:**

```json
{
  "module": "microbiome",
  "operation": "annotate",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 1,
    "text_length": 157,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "1c91805a12665e845e3dc84ecad93538c3da24850a4389154d1366268c535f9e"
  },
  "warnings": [],
  "status": "ok",
  "result": [
    {
      "id": "otu1",
      "lineage": "k__Bacteria;p__Firmicutes;c__Bacilli;o__Lactobacillales;f__Lactobacillaceae;g__Lactobacillus;s__",
      "ranks": {
        "domain": "Bacteria",
        "phylum": "Firmicutes",
        "class": "Bacilli",
        "order": "Lactobacillales",
        "family": "Lactobacillaceae",
        "genus": "Lactobacillus",
        "species": null
      },
      "lowest_rank": "genus",
      "lowest_taxon": "Lactobacillus"
    },
    {
      "id": "otu2",
      "lineage": "Archaea;Euryarchaeota;Methanobacteria",
      "ranks": {
        "domain": "Archaea",
        "phylum": "Euryarchaeota",
        "class": "Methanobacteria",
        "order": null,
        "family": null,
        "genus": null,
        "species": null
      },
      "lowest_rank": "class",
      "lowest_taxon": "Methanobacteria"
    }
  ]
}
```

## 2) assemble

- **Purpose:** Assemble metagenomic reads into contigs. Not implemented: the function always returns `status: "not_implemented"`.

- **Method:** No computation is performed and no assembler is wired into this function; a metagenome assembler such as MEGAHIT or metaSPAdes would have to be run separately. Nothing is fabricated.

- **API:** `dnakit.microbiome.assemble(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input`: any value (only summarised). `options`: ignored.

- **Output:** The call returns `status: "not_implemented"` (and `provenance.implementation` is `"not_implemented"`) with the usual envelope fields but no `result`. `required_dependency` is `a metagenome assembler (MEGAHIT, metaSPAdes) on PATH` and `message` is `Metagenome assembly is not implemented natively; no assembler is wired into this function.`

- **Example code:**

```python
import json
from dnakit.microbiome import assemble

result = assemble(["ACGTACGTACGGT", "CGTACGGTTTACA"])
print(json.dumps(result, indent=2))
```

- **Example output:**

```json
{
  "module": "microbiome",
  "operation": "assemble",
  "provenance": {
    "implementation": "not_implemented",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 2,
    "text_length": 34,
    "sequence_bases": 26,
    "lengths": [
      13,
      13
    ],
    "sha256": "b1c3a29937b1ddad6d5b46dca97b726278cd13bba542ce3c40feb88148f58dd9"
  },
  "warnings": [],
  "status": "not_implemented",
  "required_dependency": "a metagenome assembler (MEGAHIT, metaSPAdes) on PATH",
  "message": "Metagenome assembly is not implemented natively; no assembler is wired into this function."
}
```

## 3) bin

- **Purpose:** Bin metagenomic contigs into genome bins. Not implemented: the function always returns `status: "not_implemented"`.

- **Method:** No computation is performed. Binning needs read-coverage data and a binner such as MetaBAT2, MaxBin2 or CONCOCT, none of which is wired into this function. Nothing is fabricated.

- **API:** `dnakit.microbiome.bin(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input`: any value (only summarised). `options`: ignored.

- **Output:** The call returns `status: "not_implemented"` (and `provenance.implementation` is `"not_implemented"`) with the usual envelope fields but no `result`. `required_dependency` is `a contig binner (MetaBAT2, MaxBin2, CONCOCT) plus read-coverage data` and `message` is `Contig binning is not implemented natively.`

- **Example code:**

```python
import json
from dnakit.microbiome import bin as bin_contigs

result = bin_contigs({"contig1": "ACGTACGTACGGT", "contig2": "CGTACGGTTTACA"})
print(json.dumps(result, indent=2))
```

- **Example output:**

```json
{
  "module": "microbiome",
  "operation": "bin",
  "provenance": {
    "implementation": "not_implemented",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 1,
    "text_length": 56,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "d7a452d267f1cfc5abdacd676e9400fb869a1b9218a885ce2183a8ad6016f4a7"
  },
  "warnings": [],
  "status": "not_implemented",
  "required_dependency": "a contig binner (MetaBAT2, MaxBin2, CONCOCT) plus read-coverage data",
  "message": "Contig binning is not implemented natively."
}
```

## 4) build_cooccurrence_network

- **Purpose:** Build a taxon co-occurrence network from pairwise correlations of relative abundance across samples.

- **Method:** Counts are converted to per-sample relative abundance; taxa with non-zero abundance in fewer than `min_samples_present` samples are dropped. For every remaining taxon pair the correlation is Pearson's r of the relative abundances (`method="pearson"`) or of their average ranks, i.e. Spearman's rho (`"spearman"`), and a pair becomes an edge when |r| >= `min_abs_correlation`. No significance test or multiple-testing correction is applied, and correlations of relative abundances are subject to compositionality bias (the call adds a warning). At least 3 samples and at most 2000 retained taxa are required; pairs with a constant taxon (undefined r) are skipped.

- **API:** `dnakit.microbiome.build_cooccurrence_network(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input`: an abundance table: a mapping `{sample: {taxon: count}}`, `{"samples": [...], "taxa": [...], "counts": [[...]]}` (rows are samples), a list of count rows with `options["samples"]`/`options["taxa"]`, or QIIME/BIOM-TSV or CSV text/path (`options["orientation"]`: `"taxa_rows"` default | `"samples_rows"`); values must be finite and non-negative, at most 5000 samples, 200000 taxa and 20 million cells. `options`: `method` (`"spearman"` default | `"pearson"`), `min_abs_correlation` (default 0.6), `min_samples_present` (default 3), `orientation` for text tables.

- **Output:** The envelope carries `status` (`"ok"`), `module` (`"microbiome"`), `operation`, `provenance` (`implementation`, `backend`), `input_summary` (a generic text-level summary: record count, text length, base counts, SHA-256; a mapping input counts as one record), `warnings` and `result`. `result` has `method`, `min_abs_correlation`, `nodes` (retained taxa), `edges` (`source`, `target`, `correlation`, sorted by |r| descending) and `n_samples`.

- **Example code:**

```python
import json
from dnakit.microbiome import build_cooccurrence_network

table = {
    "s1": {"a": 10, "b": 20, "c": 30, "d": 1},
    "s2": {"a": 20, "b": 41, "c": 5, "d": 2},
    "s3": {"a": 30, "b": 58, "c": 25, "d": 3},
    "s4": {"a": 40, "b": 80, "c": 9, "d": 4},
}
result = build_cooccurrence_network(table, options={"method": "spearman", "min_abs_correlation": 0.9})
print(json.dumps(result, indent=2))
```

- **Example output:**

```json
{
  "module": "microbiome",
  "operation": "build_cooccurrence_network",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 1,
    "text_length": 170,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "7db91f0cb61b29d7f403f76b490c90844fc8908f32ac9d7676109e0c58527b77"
  },
  "warnings": [
    "Correlations of relative abundances are compositional; interpret with care."
  ],
  "status": "ok",
  "result": {
    "method": "spearman",
    "min_abs_correlation": 0.9,
    "nodes": [
      "a",
      "b",
      "c",
      "d"
    ],
    "edges": [
      {
        "source": "a",
        "target": "c",
        "correlation": -1.0
      },
      {
        "source": "a",
        "target": "d",
        "correlation": 1.0
      },
      {
        "source": "c",
        "target": "d",
        "correlation": -1.0
      }
    ],
    "n_samples": 4
  }
}
```

## 5) build_microbial_gene_catalog

- **Purpose:** Build a non-redundant microbial gene catalog. Not implemented: the function always returns `status: "not_implemented"`.

- **Method:** No computation is performed. A gene catalog needs gene calling (for example Prodigal) and identity clustering (CD-HIT or MMseqs2), neither of which is implemented natively. Nothing is fabricated.

- **API:** `dnakit.microbiome.build_microbial_gene_catalog(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input`: any value (only summarised). `options`: ignored.

- **Output:** The call returns `status: "not_implemented"` (and `provenance.implementation` is `"not_implemented"`) with the usual envelope fields but no `result`. `required_dependency` is `a gene predictor (Prodigal) and a sequence clusterer (CD-HIT, MMseqs2) on PATH` and `message` is `Gene-catalog construction needs gene calling and identity clustering, which are not implemented natively.`

- **Example code:**

```python
import json
from dnakit.microbiome import build_microbial_gene_catalog

result = build_microbial_gene_catalog({"contig1": "ATGAAACCCGGGTTTTAA"})
print(json.dumps(result, indent=2))
```

- **Example output:**

```json
{
  "module": "microbiome",
  "operation": "build_microbial_gene_catalog",
  "provenance": {
    "implementation": "not_implemented",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 1,
    "text_length": 33,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "546f6096d459247e3dae06ce404ad533735c07131858cf5b58f6222ad31b5654"
  },
  "warnings": [],
  "status": "not_implemented",
  "required_dependency": "a gene predictor (Prodigal) and a sequence clusterer (CD-HIT, MMseqs2) on PATH",
  "message": "Gene-catalog construction needs gene calling and identity clustering, which are not implemented natively."
}
```

## 6) build_pathogen_transmission_tree

- **Purpose:** Infer a pathogen transmission tree. Not implemented: the function always returns `status: "not_implemented"`.

- **Method:** No computation is performed. Transmission-tree inference needs a statistical (epidemiological/phylodynamic) model such as outbreaker2 or TransPhylo together with sampling dates; none is implemented. Nothing is fabricated.

- **API:** `dnakit.microbiome.build_pathogen_transmission_tree(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input`: any value (only summarised). `options`: ignored.

- **Output:** The call returns `status: "not_implemented"` (and `provenance.implementation` is `"not_implemented"`) with the usual envelope fields but no `result`. `required_dependency` is `a transmission inference model (outbreaker2, TransPhylo) with sampling dates` and `message` is `Transmission-tree inference requires a statistical model and dated samples; none is implemented.`

- **Example code:**

```python
import json
from dnakit.microbiome import build_pathogen_transmission_tree

result = build_pathogen_transmission_tree({"patientA": "ACGTACGT", "patientB": "ACGTACGA"})
print(json.dumps(result, indent=2))
```

- **Example output:**

```json
{
  "module": "microbiome",
  "operation": "build_pathogen_transmission_tree",
  "provenance": {
    "implementation": "not_implemented",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 1,
    "text_length": 48,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "e245b9769d512ad6d4af24d7600f9e72ec40ec8f60b79441a98398ef3dd81411"
  },
  "warnings": [],
  "status": "not_implemented",
  "required_dependency": "a transmission inference model (outbreaker2, TransPhylo) with sampling dates",
  "message": "Transmission-tree inference requires a statistical model and dated samples; none is implemented."
}
```

## 7) calculate_ani

- **Purpose:** Estimate average nucleotide identity between genomes from shared k-mers using the Mash formula; this is an estimate, not an alignment-based ANI.

- **Method:** For every pair of sequences the exact k-mer sets A and B are built (k-mers containing non-ACGT symbols are skipped; with `canonical` a k-mer and its reverse complement are one word). From the Jaccard index j = |A&B| / |A|B| the Mash distance is D = -ln(2j / (1 + j)) / k and `ani_estimate` = 1 - D (`null` when no k-mer is shared); the containments |A&B|/|A| and |A&B|/|B| are also reported. Unlike Mash itself no MinHash sketch is used; the full k-mer sets are compared. This is not FastANI or ANIb (there is no fragment alignment), and the estimate is only reliable for genomes of similar size with ANI above about 0.8; it is not validated against alignment-based ANI here. The call adds a warning to this effect.

- **API:** `dnakit.microbiome.calculate_ani(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input`: DNA sequences as a list of strings, a mapping `{id: sequence}` or FASTA text/path (ACGT plus IUPAC codes; U is read as T); between 2 and 50 sequences, at most 20 Mb in total. `options`: `k` (default 21, 2-32), `canonical` (default `true`).

- **Output:** The envelope carries `status` (`"ok"`), `module` (`"microbiome"`), `operation`, `provenance` (`implementation`, `backend`), `input_summary` (a generic text-level summary: record count, text length, base counts, SHA-256; a mapping input counts as one record), `warnings` and `result`. `result` has `method` (`"mash_kmer_jaccard"`) and `pairs` (for each unordered pair `a`, `b`, `k`, `shared_kmers`, `jaccard`, `containment_a_in_b`, `containment_b_in_a`, `mash_distance`, `ani_estimate`).

- **Example code:**

```python
import json
import random
from dnakit.microbiome import calculate_ani

rng = random.Random(1)
genome_a = "".join(rng.choice("ACGT") for _ in range(3000))
genome_b = "".join(rng.choice("ACGT".replace(c, "")) if i % 50 == 0 else c   # change every 50th base
                   for i, c in enumerate(genome_a))
result = calculate_ani({"A": genome_a, "B": genome_b}, options={"k": 21})
print(json.dumps(result, indent=2))
```

- **Example output:**

```json
{
  "module": "microbiome",
  "operation": "calculate_ani",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 1,
    "text_length": 6018,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "b5509626509603ca23fb0e0c995cf060a54a720b39e6e6a3a06732cf2ab35eaf"
  },
  "warnings": [
    "k-mer set (Mash-formula) ANI estimate, not an alignment-based ANI."
  ],
  "status": "ok",
  "result": {
    "method": "mash_kmer_jaccard",
    "pairs": [
      {
        "a": "A",
        "b": "B",
        "k": 21,
        "shared_kmers": 1740.0,
        "jaccard": 0.41232227488151657,
        "containment_a_in_b": 0.5838926174496645,
        "containment_b_in_a": 0.5838926174496645,
        "mash_distance": 0.025620866061470255,
        "ani_estimate": 0.9743791339385297
      }
    ]
  }
}
```

## 8) calculate_differential_abundance

- **Purpose:** Test every taxon for differential relative abundance between two sample groups (Mann-Whitney U, Benjamini-Hochberg adjusted).

- **Method:** Counts are converted to per-sample relative abundance and taxa present in a fraction of samples below `min_prevalence` are skipped. For each remaining taxon the group means are computed, `log2_fold_change` = log2((mean_b + pseudocount) / (mean_a + pseudocount)), and the two-sided Mann-Whitney U test of group A against group B is run with `scipy.stats.mannwhitneyu` (U is reported for group A; if all values are identical U = n_a*n_b/2 and p = 1). The p-values of all tested taxa are Benjamini-Hochberg adjusted to `q_value`, and taxa are sorted by p-value. This is a rank test on compositional relative abundances: no normalisation beyond relative abundance, no covariates and no zero-inflation model. scipy is required: without it the function returns `status: "not_implemented"` with `required_dependency` `pip install scipy`.

- **API:** `dnakit.microbiome.calculate_differential_abundance(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input`: an abundance table: a mapping `{sample: {taxon: count}}`, `{"samples": [...], "taxa": [...], "counts": [[...]]}` (rows are samples), a list of count rows with `options["samples"]`/`options["taxa"]`, or QIIME/BIOM-TSV or CSV text/path (`options["orientation"]`: `"taxa_rows"` default | `"samples_rows"`); values must be finite and non-negative, at most 5000 samples, 200000 taxa and 20 million cells. `options`: `groups` (required: a list of group labels in sample order, or a mapping sample -> label); `group_a`/`group_b` (default: the two sorted labels, which must then be exactly two); `pseudocount` (default 1e-6); `min_prevalence` (default 0.1). Each group needs at least two samples.

- **Output:** The envelope carries `status` (`"ok"`), `module` (`"microbiome"`), `operation`, `provenance` (`implementation`, `backend`), `input_summary` (a generic text-level summary: record count, text length, base counts, SHA-256; a mapping input counts as one record), `warnings` and `result`. `result` has `group_a`, `group_b`, `n_a`, `n_b`, `test` (`"mann_whitney_u"`) and `taxa` (per taxon `taxon`, `mean_a`, `mean_b`, `log2_fold_change`, `u_statistic`, `p_value`, `q_value`).

- **Example code:**

```python
import json
from dnakit.microbiome import calculate_differential_abundance

table = {
    "c1": {"a": 90, "b": 10, "c": 50}, "c2": {"a": 80, "b": 20, "c": 55},
    "c3": {"a": 85, "b": 15, "c": 45},
    "t1": {"a": 20, "b": 80, "c": 50}, "t2": {"a": 30, "b": 70, "c": 52},
    "t3": {"a": 25, "b": 75, "c": 48},
}
groups = {"c1": "control", "c2": "control", "c3": "control", "t1": "treated", "t2": "treated", "t3": "treated"}
result = calculate_differential_abundance(table, options={"groups": groups})
print(json.dumps(result, indent=2))
```

- **Example output:**

```json
{
  "module": "microbiome",
  "operation": "calculate_differential_abundance",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 1,
    "text_length": 210,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "3d66b69cc2ad0a8cda972aa64d7d905e73c3c44f81aff7edb99364aaa3b4ad65"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "group_a": "control",
    "group_b": "treated",
    "n_a": 3,
    "n_b": 3,
    "taxa": [
      {
        "taxon": "a",
        "mean_a": 0.5674453096032629,
        "mean_b": 0.16654022443496128,
        "log2_fold_change": -1.7686045585338197,
        "u_statistic": 9.0,
        "p_value": 0.1,
        "q_value": 0.15000000000000002
      },
      {
        "taxon": "b",
        "mean_a": 0.09971573353108391,
        "mean_b": 0.5002054686265213,
        "log2_fold_change": 2.326616184784275,
        "u_statistic": 0.0,
        "p_value": 0.1,
        "q_value": 0.15000000000000002
      },
      {
        "taxon": "c",
        "mean_a": 0.3328389568656532,
        "mean_b": 0.3332543069385175,
        "log2_fold_change": 0.0017992128170873996,
        "u_statistic": 4.5,
        "p_value": 1.0,
        "q_value": 1.0
      }
    ],
    "test": "mann_whitney_u"
  }
}
```

## 9) calculate_digital_ddh

- **Purpose:** Compute digital DNA-DNA hybridisation (dDDH). Not implemented: the function always returns `status: "not_implemented"`.

- **Method:** No computation is performed. dDDH needs genome alignments (BLAST, BLAT or MUMmer HSPs) and the GGDC regression formulas, which are not implemented. For a k-mer based identity estimate use `calculate_ani`.

- **API:** `dnakit.microbiome.calculate_digital_ddh(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input`: any value (only summarised). `options`: ignored.

- **Output:** The call returns `status: "not_implemented"` (and `provenance.implementation` is `"not_implemented"`) with the usual envelope fields but no `result`. `required_dependency` is `genome alignments (BLAST/BLAT/MUMmer HSPs) and the GGDC regression formulas` and `message` is `dDDH is not implemented; use calculate_ani for a k-mer identity estimate instead.`

- **Example code:**

```python
import json
from dnakit.microbiome import calculate_digital_ddh

result = calculate_digital_ddh({"A": "ACGTACGTACGT", "B": "ACGTACGTACGA"})
print(json.dumps(result, indent=2))
```

- **Example output:**

```json
{
  "module": "microbiome",
  "operation": "calculate_digital_ddh",
  "provenance": {
    "implementation": "not_implemented",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 1,
    "text_length": 42,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "18186543adfa540845c61ad49fbb8968a47e742c04b3030c93180cc48804768a"
  },
  "warnings": [],
  "status": "not_implemented",
  "required_dependency": "genome alignments (BLAST/BLAT/MUMmer HSPs) and the GGDC regression formulas",
  "message": "dDDH is not implemented; use calculate_ani for a k-mer identity estimate instead."
}
```

## 10) classify

- **Purpose:** Assign taxonomy to reads or contigs. Not implemented natively: the function always returns `status: "not_implemented"`.

- **Method:** No computation is performed because there is no reference database. The message points to the `classify_*` wrappers (for example `classify_kraken2`, `classify_kaiju`, `classify_centrifuge`, `classify_gtdbtk`) of the same module, which call an installed classifier and database. Nothing is fabricated.

- **API:** `dnakit.microbiome.classify(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input`: any value (only summarised). `options`: ignored.

- **Output:** The call returns `status: "not_implemented"` (and `provenance.implementation` is `"not_implemented"`) with the usual envelope fields but no `result`. `required_dependency` is `a taxonomic classifier database (Kraken2, Kaiju, Centrifuge, GTDB-Tk; see classify_kraken2 etc.)` and `message` is `Native taxonomic classification is not implemented; use the classify_* wrappers with an installed database.`

- **Example code:**

```python
import json
from dnakit.microbiome import classify

result = classify(["ACGTACGTACGTACGT"])
print(json.dumps(result, indent=2))
```

- **Example output:**

```json
{
  "module": "microbiome",
  "operation": "classify",
  "provenance": {
    "implementation": "not_implemented",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 1,
    "text_length": 20,
    "sequence_bases": 16,
    "lengths": [
      16
    ],
    "sha256": "87cd4b17356c3dee4b13d54fba3beefea5678a7b651616a23b27c782da631af2"
  },
  "warnings": [],
  "status": "not_implemented",
  "required_dependency": "a taxonomic classifier database (Kraken2, Kaiju, Centrifuge, GTDB-Tk; see classify_kraken2 etc.)",
  "message": "Native taxonomic classification is not implemented; use the classify_* wrappers with an installed database."
}
```

## 11) compare

- **Purpose:** Compute a between-sample beta-diversity distance matrix from an abundance table.

- **Method:** For every pair of samples a, b: `braycurtis` = sum|a-b| / sum(a+b) on the (raw or relative) counts; `jaccard` = 1 - shared/union of presence (count > 0); `aitchison` = Euclidean distance between centred log-ratio profiles of counts + `pseudocount`; `euclidean` = Euclidean distance of the values. With `relative` every sample is first converted to relative abundance. The matrix is symmetric with a zero diagonal; no rarefaction or phylogenetic metric (UniFrac) is used.

- **API:** `dnakit.microbiome.compare(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input`: an abundance table: a mapping `{sample: {taxon: count}}`, `{"samples": [...], "taxa": [...], "counts": [[...]]}` (rows are samples), a list of count rows with `options["samples"]`/`options["taxa"]`, or QIIME/BIOM-TSV or CSV text/path (`options["orientation"]`: `"taxa_rows"` default | `"samples_rows"`); values must be finite and non-negative, at most 5000 samples, 200000 taxa and 20 million cells; at least 2 and at most 2000 samples. `options`: `entity` (only `"beta_diversity"`); `metric` (`"braycurtis"` default | `"jaccard"` | `"aitchison"` | `"euclidean"`); `pseudocount` (default 1.0, used by `aitchison`); `relative` (default `false`).

- **Output:** The envelope carries `status` (`"ok"`), `module` (`"microbiome"`), `operation`, `provenance` (`implementation`, `backend`), `input_summary` (a generic text-level summary: record count, text length, base counts, SHA-256; a mapping input counts as one record), `warnings` and `result`. `result` has `metric`, `samples` and the symmetric `distance_matrix` (rows and columns follow `samples`).

- **Example code:**

```python
import json
from dnakit.microbiome import compare

table = {"s1": {"a": 5, "b": 5, "c": 5, "d": 5},
         "s2": {"a": 1, "b": 1, "c": 2, "d": 6},
         "s3": {"a": 10, "b": 0, "c": 3, "d": 0}}
result = compare(table, options={"metric": "braycurtis"})
print(json.dumps(result, indent=2))
```

- **Example output:**

```json
{
  "module": "microbiome",
  "operation": "compare",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 1,
    "text_length": 121,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "4193a3bec6078b3f816617bfa2ee0d04d6170d8243e9c6557184867a5162df81"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "metric": "braycurtis",
    "samples": [
      "s1",
      "s2",
      "s3"
    ],
    "distance_matrix": [
      [
        0.0,
        0.4,
        0.5151515151515151
      ],
      [
        0.4,
        0.0,
        0.7391304347826086
      ],
      [
        0.5151515151515151,
        0.7391304347826086,
        0.0
      ]
    ]
  }
}
```

## 12) dereplicate

- **Purpose:** Collapse redundant genome or sequence records into clusters with one representative each.

- **Method:** `method="exact"` groups identical sequences (case-insensitive; with `canonical` also sequences identical up to reverse complement). `method="ani"` is greedy: sequences are sorted longest first (then by id), and each one joins the first existing cluster whose representative has a k-mer ANI estimate >= `min_ani` with it, otherwise it starts a new cluster. The ANI estimate is the Mash-formula estimate of `calculate_ani` computed on exact k-mer sets (`k` default 16 here), so the `"ani"` mode is a k-mer estimate and not an alignment-based clustering such as dRep; it also compares only with each cluster's representative.

- **API:** `dnakit.microbiome.dereplicate(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input`: DNA sequences as a list of strings, a mapping `{id: sequence}` or FASTA text/path (ACGT plus IUPAC codes; U is read as T); 1 to 2000 sequences, at most 50 Mb in total. `options`: `method` (`"exact"` default | `"ani"`), `canonical` (default `true`), `k` (default 16, 2-32, `ani` mode), `min_ani` (default 0.95, `ani` mode).

- **Output:** The envelope carries `status` (`"ok"`), `module` (`"microbiome"`), `operation`, `provenance` (`implementation`, `backend`), `input_summary` (a generic text-level summary: record count, text length, base counts, SHA-256; a mapping input counts as one record), `warnings` and `result`. `result` has `method`, `input_count`, `cluster_count` and `clusters` (each with `representative` and `members`).

- **Example code:**

```python
import json
import random
from dnakit.microbiome import dereplicate

rng = random.Random(2)
g1 = "".join(rng.choice("ACGT") for _ in range(2000))
g2 = "".join(rng.choice("ACGT".replace(c, "")) if i % 200 == 0 else c for i, c in enumerate(g1))   # ~0.5 % changed
g3 = "".join(rng.choice("ACGT") for _ in range(1800))                                              # unrelated
result = dereplicate({"g1": g1, "g2": g2, "g3": g3}, options={"method": "ani", "k": 16, "min_ani": 0.95})
print(json.dumps(result, indent=2))
```

- **Example output:**

```json
{
  "module": "microbiome",
  "operation": "dereplicate",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 1,
    "text_length": 5830,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "340d4746b48c4e13c53b16f2046187b727af0cd4ee1b052ad4cea41fe5176b8e"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "method": "ani",
    "input_count": 3,
    "cluster_count": 2,
    "clusters": [
      {
        "representative": "g1",
        "members": [
          "g1",
          "g2"
        ]
      },
      {
        "representative": "g3",
        "members": [
          "g3"
        ]
      }
    ]
  }
}
```

## 13) detect

- **Purpose:** Detect taxon prevalence and the core microbiome (taxa present above a detection level in enough samples).

- **Method:** Relative abundance of a taxon in a sample is count / sample total. A taxon counts as present in a sample when its relative abundance is strictly greater than `detection`; prevalence is the fraction of samples where it is present. Core taxa have prevalence >= `prevalence`. Taxa are sorted by prevalence, then mean relative abundance. No statistical test or abundance-weighted core definition is used.

- **API:** `dnakit.microbiome.detect(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input`: an abundance table: a mapping `{sample: {taxon: count}}`, `{"samples": [...], "taxa": [...], "counts": [[...]]}` (rows are samples), a list of count rows with `options["samples"]`/`options["taxa"]`, or QIIME/BIOM-TSV or CSV text/path (`options["orientation"]`: `"taxa_rows"` default | `"samples_rows"`); values must be finite and non-negative, at most 5000 samples, 200000 taxa and 20 million cells; every sample needs a positive total. `options`: `entity` (`"core_microbiome"` default, adds the core taxa | `"prevalence"`, only the table); `detection` (default 0.0); `prevalence` (default 0.5).

- **Output:** The envelope carries `status` (`"ok"`), `module` (`"microbiome"`), `operation`, `provenance` (`implementation`, `backend`), `input_summary` (a generic text-level summary: record count, text length, base counts, SHA-256; a mapping input counts as one record), `warnings` and `result`. `result` has `detection`, `n_samples` and `prevalence_table` (per taxon `taxon`, `samples_present`, `prevalence`, `mean_relative_abundance`); for `core_microbiome` also `prevalence_threshold` and `core_taxa`.

- **Example code:**

```python
import json
from dnakit.microbiome import detect

table = {"s1": {"a": 50, "b": 10, "c": 0}, "s2": {"a": 40, "b": 0, "c": 0}, "s3": {"a": 30, "b": 5, "c": 1}}
result = detect(table, options={"entity": "core_microbiome", "detection": 0.0, "prevalence": 0.6})
print(json.dumps(result, indent=2))
```

- **Example output:**

```json
{
  "module": "microbiome",
  "operation": "detect",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 1,
    "text_length": 100,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "417c10ac425501d532651729f25cccaad64f71d0508e4209810a058e285d982c"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "detection": 0.0,
    "n_samples": 3,
    "prevalence_table": [
      {
        "taxon": "a",
        "samples_present": 3,
        "prevalence": 1.0,
        "mean_relative_abundance": 0.888888888888889
      },
      {
        "taxon": "b",
        "samples_present": 2,
        "prevalence": 0.6666666666666666,
        "mean_relative_abundance": 0.10185185185185186
      },
      {
        "taxon": "c",
        "samples_present": 1,
        "prevalence": 0.3333333333333333,
        "mean_relative_abundance": 0.009259259259259259
      }
    ],
    "prevalence_threshold": 0.6,
    "core_taxa": [
      "a",
      "b"
    ]
  }
}
```

## 14) export

- **Purpose:** Serialise an abundance table as QIIME-style TSV (taxa as rows) or JSON, optionally writing it to a file.

- **Method:** The table is parsed and validated, then written as TSV with the header `#OTU ID` followed by the sample ids and one row per taxon (whole numbers are written without a decimal point), or as JSON of `samples`, `taxa` and `counts` with sorted keys. An existing output file is replaced only when `overwrite` is true. The text is also returned in `result`.

- **API:** `dnakit.microbiome.export(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input`: an abundance table: a mapping `{sample: {taxon: count}}`, `{"samples": [...], "taxa": [...], "counts": [[...]]}` (rows are samples), a list of count rows with `options["samples"]`/`options["taxa"]`, or QIIME/BIOM-TSV or CSV text/path (`options["orientation"]`: `"taxa_rows"` default | `"samples_rows"`); values must be finite and non-negative, at most 5000 samples, 200000 taxa and 20 million cells. `options`: `format` (`"tsv"` default | `"json"`); `output` (optional file path); `overwrite` (default `false`).

- **Output:** The envelope carries `status` (`"ok"`), `module` (`"microbiome"`), `operation`, `provenance` (`implementation`, `backend`), `input_summary` (a generic text-level summary: record count, text length, base counts, SHA-256; a mapping input counts as one record), `warnings` and `result`. `result` has `format`, `text` and `output_path` (the written file or `null`).

- **Example code:**

```python
import json
from dnakit.microbiome import export

table = {"s1": {"a": 5, "b": 0}, "s2": {"a": 1, "b": 7}}
result = export(table, options={"format": "tsv"})
print(json.dumps(result, indent=2))
```

- **Example output:**

```json
{
  "module": "microbiome",
  "operation": "export",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 1,
    "text_length": 48,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "826172e6734fde444cb3c4c03ccd74ba350458e29235dc64e5827464038ab612"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "format": "tsv",
    "text": "#OTU ID\ts1\ts2\na\t5\t1\nb\t0\t7\n",
    "output_path": null
  }
}
```

## 15) infer_viral_quasispecies

- **Purpose:** Reconstruct viral quasispecies haplotypes. Not implemented: the function always returns `status: "not_implemented"`.

- **Method:** No computation is performed. Haplotype reconstruction needs aligned reads and a program such as QuRe, ShoRAH or PredictHaplo; none is wired into this function. Nothing is fabricated.

- **API:** `dnakit.microbiome.infer_viral_quasispecies(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input`: any value (only summarised). `options`: ignored.

- **Output:** The call returns `status: "not_implemented"` (and `provenance.implementation` is `"not_implemented"`) with the usual envelope fields but no `result`. `required_dependency` is `a viral haplotype reconstruction tool (QuRe, ShoRAH, PredictHaplo) and aligned reads` and `message` is `Quasispecies inference is not implemented natively.`

- **Example code:**

```python
import json
from dnakit.microbiome import infer_viral_quasispecies

result = infer_viral_quasispecies(["ACGTACGTAC", "ACGTACGTAA"])
print(json.dumps(result, indent=2))
```

- **Example output:**

```json
{
  "module": "microbiome",
  "operation": "infer_viral_quasispecies",
  "provenance": {
    "implementation": "not_implemented",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 2,
    "text_length": 28,
    "sequence_bases": 20,
    "lengths": [
      10,
      10
    ],
    "sha256": "87931cb5c372cda553027980d7b4eb79ea0a83bd1432f85ad02bba23469f8907"
  },
  "warnings": [],
  "status": "not_implemented",
  "required_dependency": "a viral haplotype reconstruction tool (QuRe, ShoRAH, PredictHaplo) and aligned reads",
  "message": "Quasispecies inference is not implemented natively."
}
```

## 16) plot

- **Purpose:** Draw a stacked relative-abundance bar chart (one bar per sample) as an SVG string.

- **Method:** Relative abundance per sample is drawn as stacked segments in a hand-written SVG (no plotting library). The `top_n` taxa with the largest total count over all samples get their own colour from a fixed 10-colour palette and all other taxa are pooled into a grey `Other` segment; a legend lists the top taxa. Each segment is a `<rect class="bar">` carrying `data-sample`, `data-taxon` and `data-fraction` attributes.

- **API:** `dnakit.microbiome.plot(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input`: an abundance table: a mapping `{sample: {taxon: count}}`, `{"samples": [...], "taxa": [...], "counts": [[...]]}` (rows are samples), a list of count rows with `options["samples"]`/`options["taxa"]`, or QIIME/BIOM-TSV or CSV text/path (`options["orientation"]`: `"taxa_rows"` default | `"samples_rows"`); values must be finite and non-negative, at most 5000 samples, 200000 taxa and 20 million cells; every sample needs a positive total. `options`: `top_n` (default 8, range 1-10); `output` (optional `.svg` path).

- **Output:** The envelope carries `status` (`"ok"`), `module` (`"microbiome"`), `operation`, `provenance` (`implementation`, `backend`), `input_summary` (a generic text-level summary: record count, text length, base counts, SHA-256; a mapping input counts as one record), `warnings` and `result`. `result` has `svg` (the SVG text) and `output_path` (the written file or `null`).

- **Example code:**

```python
import json
from dnakit.microbiome import plot

table = {"s1": {"a": 8, "b": 2}, "s2": {"a": 3, "b": 7}}
result = plot(table, options={"top_n": 2})
print(json.dumps(result, indent=2))
```

- **Example output:** Long arrays and strings in the output are shortened with `...`.

```json
{
  "module": "microbiome",
  "operation": "plot",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 1,
    "text_length": 48,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "ff58a93da8baaf87e71302a304526230cd06541032b0298e374387dc98fd9d4b"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "svg": "<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"312\" height=\"280\"><rect class=\"bar\" data-sample=\"s1\" data-taxon=\"a\" data-fraction=\"...",
    "output_path": null
  }
}
```

## 17) predict

- **Purpose:** Predict microbiome phenotypes or functions. Not implemented: the function always returns `status: "not_implemented"`.

- **Method:** No computation is performed. Predicting function or phenotype from a profile needs trained reference data such as PICRUSt2 or BugBase, which is not bundled. Nothing is fabricated.

- **API:** `dnakit.microbiome.predict(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input`: any value (only summarised). `options`: ignored.

- **Output:** The call returns `status: "not_implemented"` (and `provenance.implementation` is `"not_implemented"`) with the usual envelope fields but no `result`. `required_dependency` is `a trained model or reference database (e.g. PICRUSt2, BugBase)` and `message` is `Prediction from a microbiome profile needs trained reference data that is not bundled.`

- **Example code:**

```python
import json
from dnakit.microbiome import predict

result = predict({"s1": {"a": 8, "b": 2}})
print(json.dumps(result, indent=2))
```

- **Example output:**

```json
{
  "module": "microbiome",
  "operation": "predict",
  "provenance": {
    "implementation": "not_implemented",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 1,
    "text_length": 24,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "beafe4e6da1b369b3fbf3192253bb6e2e35f1331de917590e4047cfcb11ee9ce"
  },
  "warnings": [],
  "status": "not_implemented",
  "required_dependency": "a trained model or reference database (e.g. PICRUSt2, BugBase)",
  "message": "Prediction from a microbiome profile needs trained reference data that is not bundled."
}
```

## 18) qc

- **Purpose:** Quality-check an abundance table: per-sample depth, richness and zero fraction plus low-depth flags.

- **Method:** For each sample: `total` (sum of counts), `observed` (taxa with count > 0), `zero_fraction` = 1 - observed / number of taxa, and `low_depth` when the total is below `min_depth`. Taxa that are zero in every sample are listed as `empty_taxa`. The depth threshold is a user-chosen rule of thumb, not derived from the data.

- **API:** `dnakit.microbiome.qc(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input`: an abundance table: a mapping `{sample: {taxon: count}}`, `{"samples": [...], "taxa": [...], "counts": [[...]]}` (rows are samples), a list of count rows with `options["samples"]`/`options["taxa"]`, or QIIME/BIOM-TSV or CSV text/path (`options["orientation"]`: `"taxa_rows"` default | `"samples_rows"`); values must be finite and non-negative, at most 5000 samples, 200000 taxa and 20 million cells. `options`: `entity` (only `"abundance_table"`); `min_depth` (default 1000).

- **Output:** The envelope carries `status` (`"ok"`), `module` (`"microbiome"`), `operation`, `provenance` (`implementation`, `backend`), `input_summary` (a generic text-level summary: record count, text length, base counts, SHA-256; a mapping input counts as one record), `warnings` and `result`. `result` has `n_samples`, `n_taxa`, `min_depth`, `samples` (per sample `sample`, `total`, `observed`, `zero_fraction`, `low_depth`), `empty_taxa` and `flagged_samples`.

- **Example code:**

```python
import json
from dnakit.microbiome import qc

table = {"s1": {"a": 900, "b": 100, "c": 0}, "s2": {"a": 30, "b": 0, "c": 0}}
result = qc(table, options={"min_depth": 1000})
print(json.dumps(result, indent=2))
```

- **Example output:**

```json
{
  "module": "microbiome",
  "operation": "qc",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 1,
    "text_length": 69,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "3f254bebd28a310387a75982cd4bf7766e8f0eb10be64aee81359716e2eb18ba"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "n_samples": 2,
    "n_taxa": 3,
    "min_depth": 1000.0,
    "samples": [
      {
        "sample": "s1",
        "total": 1000.0,
        "observed": 2,
        "zero_fraction": 0.33333333333333337,
        "low_depth": false
      },
      {
        "sample": "s2",
        "total": 30.0,
        "observed": 1,
        "zero_fraction": 0.6666666666666667,
        "low_depth": true
      }
    ],
    "empty_taxa": [
      "c"
    ],
    "flagged_samples": [
      "s2"
    ]
  }
}
```

## 19) quantify

- **Purpose:** Quantify an abundance table: relative abundance, alpha diversity, rarefaction or centred log-ratio (CLR). Reads are not mapped here; counts must already exist.

- **Method:** `relative_abundance` (alias `abundance`): count / sample total. `alpha_diversity`, per sample: observed taxa, Shannon H = -sum p ln p (natural log unless `log_base`), Simpson (Gini-Simpson) = 1 - sum p^2, inverse Simpson = 1 / sum p^2, Pielou evenness = H / ln(observed) (`null` for one taxon), bias-corrected Chao1 = S_obs + F1 (F1 - 1) / (2 (F2 + 1)) and Good's coverage = 1 - F1 / N, where F1/F2 are the numbers of singleton/doubleton taxa and N the total count (Chao1 and coverage only for integer counts, otherwise `null`). `rarefaction`: each sample is subsampled to `depth` reads without replacement with `random.Random(seed)`; samples with fewer reads are dropped and listed (or raise an error with `drop_short` false). `clr`: log(count + `pseudocount`) minus its mean over taxa.

- **API:** `dnakit.microbiome.quantify(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input`: an abundance table: a mapping `{sample: {taxon: count}}`, `{"samples": [...], "taxa": [...], "counts": [[...]]}` (rows are samples), a list of count rows with `options["samples"]`/`options["taxa"]`, or QIIME/BIOM-TSV or CSV text/path (`options["orientation"]`: `"taxa_rows"` default | `"samples_rows"`); values must be finite and non-negative, at most 5000 samples, 200000 taxa and 20 million cells; every sample needs a positive total except for `rarefaction`. `options`: `entity` (`"relative_abundance"` default | `"abundance"` | `"alpha_diversity"` | `"rarefaction"` | `"clr"`); `log_base` (default e, `alpha_diversity`); `depth` (required for `rarefaction`, integer counts only), `seed` (default 0), `drop_short` (default `true`); `pseudocount` (default 1.0, `clr`).

- **Output:** The envelope carries `status` (`"ok"`), `module` (`"microbiome"`), `operation`, `provenance` (`implementation`, `backend`), `input_summary` (a generic text-level summary: record count, text length, base counts, SHA-256; a mapping input counts as one record), `warnings` and `result`. `result` depends on `entity`: `relative_abundance` gives `samples`, `taxa`, `relative_abundance`; `alpha_diversity` gives `samples` and `metrics` (per sample `observed`, `shannon`, `simpson`, `inverse_simpson`, `pielou_evenness`, `chao1`, `goods_coverage`); `rarefaction` gives `depth`, `dropped_samples`, `samples`, `taxa`, `counts`; `clr` gives `samples`, `taxa`, `clr`.

- **Example code:**

```python
import json
from dnakit.microbiome import quantify

table = {"s1": {"a": 5, "b": 5, "c": 5, "d": 5}, "s2": {"a": 1, "b": 1, "c": 2, "d": 6}}
result = quantify(table, options={"entity": "alpha_diversity"})
print(json.dumps(result, indent=2))
```

- **Example output:**

```json
{
  "module": "microbiome",
  "operation": "quantify",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 1,
    "text_length": 80,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "ea761b00be122b589f5902ee4cb3fb88280b57aad31592d0352b9e8c8102d34f"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "samples": [
      "s1",
      "s2"
    ],
    "metrics": {
      "s1": {
        "observed": 4.0,
        "shannon": 1.3862943611198906,
        "simpson": 0.75,
        "inverse_simpson": 4.0,
        "pielou_evenness": 1.0,
        "chao1": 4.0,
        "goods_coverage": 1.0
      },
      "s2": {
        "observed": 4.0,
        "shannon": 1.0888999753452238,
        "simpson": 0.5800000000000001,
        "inverse_simpson": 2.380952380952381,
        "pielou_evenness": 0.7854752972273344,
        "chao1": 4.5,
        "goods_coverage": 0.8
      }
    }
  }
}
```

## 20) read

- **Purpose:** Read an abundance table (QIIME/BIOM-TSV or CSV text, a file path or an in-memory structure) into a samples x taxa count matrix.

- **Method:** Text is split at tabs (or commas); the first line that is not a `#` comment (a header starting `#OTU`, `#Feature` or `#ID` counts as the header) gives the sample ids, and the first column gives the taxon ids. With `taxa_rows` the table is transposed so that rows of the result are samples. Values must be finite and non-negative numbers, ids must be unique and all rows must have the same width; otherwise an input error is raised. No normalisation is applied.

- **API:** `dnakit.microbiome.read(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input`: an abundance table: a mapping `{sample: {taxon: count}}`, `{"samples": [...], "taxa": [...], "counts": [[...]]}` (rows are samples), a list of count rows with `options["samples"]`/`options["taxa"]`, or QIIME/BIOM-TSV or CSV text/path (`options["orientation"]`: `"taxa_rows"` default | `"samples_rows"`); values must be finite and non-negative, at most 5000 samples, 200000 taxa and 20 million cells. `options`: `orientation` for text tables.

- **Output:** The envelope carries `status` (`"ok"`), `module` (`"microbiome"`), `operation`, `provenance` (`implementation`, `backend`), `input_summary` (a generic text-level summary: record count, text length, base counts, SHA-256; a mapping input counts as one record), `warnings` and `result`. `result` has `samples`, `taxa` and `counts` (rows are samples, values are floats).

- **Example code:**

```python
import json
from dnakit.microbiome import read

text = "#OTU ID\ts1\ts2\nEscherichia\t10\t0\nBacteroides\t5\t20\n"
result = read(text)
print(json.dumps(result, indent=2))
```

- **Example output:**

```json
{
  "module": "microbiome",
  "operation": "read",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 3,
    "text_length": 48,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "a3f580da762132a4f472fe7656ecdcf87576a39c444c4a62471110133decc416"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "samples": [
      "s1",
      "s2"
    ],
    "taxa": [
      "Escherichia",
      "Bacteroides"
    ],
    "counts": [
      [
        10.0,
        5.0
      ],
      [
        0.0,
        20.0
      ]
    ]
  }
}
```

## 21) track

- **Purpose:** Track strains or pathogens across samples. Not implemented: the function always returns `status: "not_implemented"`.

- **Method:** No computation is performed. Strain tracking needs strain-level marker or SNP references (for example StrainPhlAn, inStrain) and epidemiological metadata, which are not available natively. Nothing is fabricated.

- **API:** `dnakit.microbiome.track(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input`: any value (only summarised). `options`: ignored.

- **Output:** The call returns `status: "not_implemented"` (and `provenance.implementation` is `"not_implemented"`) with the usual envelope fields but no `result`. `required_dependency` is `strain-level marker/SNP references (e.g. StrainPhlAn, inStrain) and epidemiological metadata` and `message` is `Strain/pathogen tracking is not implemented natively.`

- **Example code:**

```python
import json
from dnakit.microbiome import track

result = track({"isolate1": "ACGTACGT", "isolate2": "ACGTACGA"})
print(json.dumps(result, indent=2))
```

- **Example output:**

```json
{
  "module": "microbiome",
  "operation": "track",
  "provenance": {
    "implementation": "not_implemented",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 1,
    "text_length": 48,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "130efd49bcf3ffab1415417f4a23967523f5c89a74650caf0c6f17f20e766ba1"
  },
  "warnings": [],
  "status": "not_implemented",
  "required_dependency": "strain-level marker/SNP references (e.g. StrainPhlAn, inStrain) and epidemiological metadata",
  "message": "Strain/pathogen tracking is not implemented natively."
}
```
