# Haplotypes

This category contains **17** documented capabilities; 16 compute a result natively in pure Python from the supplied genotypes or haplotypes (trio and pedigree logic, haplotype statistics, LD blocks and tag SNPs, EM haplotype reconstruction, simulation, SVG); `count` is a generic record counter that counts non-blank lines or list items and does not interpret them as haplotypes. Only `impute` always returns `status: "not_implemented"`, because it needs a phased reference panel and an imputation engine that are not part of the input. No external program or optional library is required. The documented functions are exported through the formal DNAKit namespace; functions with a compatible JSON call plan are also exposed through the Agent/MCP registry.

## 1) analyze_trio_inheritance

- **Purpose:** Checks the Mendelian consistency of a child's genotypes against both parents and reports the transmitted alleles.

- **Method:** At each site the child's two alleles must be formable as one allele from the father plus one allele from the mother; otherwise the site is `inconsistent` (a Mendelian error). A site with any missing call is `missing` and skipped. When exactly one (paternal, maternal) allele assignment explains the child's genotype it is reported; when several do (for example all three are heterozygous) the origin is left undetermined. Mendelian error rate = inconsistent sites / informative sites (consistent + inconsistent). Autosomal, diploid genotypes are assumed; no genotyping-error model and no X-chromosome handling.

- **API:** `dnakit.haplotypes.analyze_trio_inheritance(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input` (required): `{"child": [...], "father": [...], "mother": [...]}` with one genotype call per site (lists of equal length) and an optional `"sites"` list of labels. A call is `"a/b"` or `"a|b"` with integer alleles, a pair `(a, b)`, an alt-allele dosage 0/1/2, or `None` / `"./."` for missing. `options`: none.

- **Output:** Standard envelope: `status` (`"ok"`), `module`, `operation`, `provenance`, `input_summary`, `warnings` and `result`. `result` contains `site_count`, `informative_sites`, `consistent`, `inconsistent`, `missing`, `mendelian_error_rate` (null if no informative site), `inconsistent_indices`, `sites` (one entry per site with `index`, optional `site` label, `status` and, when determinable, `paternal_allele` / `maternal_allele`) and a `rule` string.

- **Example code:**

```python
import json

from dnakit.haplotypes import analyze_trio_inheritance

trio = {
    "sites": ["rs1", "rs2", "rs3", "rs4"],
    "child": ["0/1", "1/1", "0/0", "./."],
    "father": ["0/0", "0/1", "0/1", "0/0"],
    "mother": ["0/1", "0/1", "1/1", "0/0"],
}
result = analyze_trio_inheritance(trio)
print(json.dumps(result, indent=2))
```

- **Example output:**

```json
{
  "module": "haplotypes",
  "operation": "analyze_trio_inheritance",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 1,
    "text_length": 158,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "dec8cb6e39f7932d76a0417af1afaff5105efb3538d6d64a2e3718957d0cb659"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "site_count": 4,
    "informative_sites": 3,
    "consistent": 2,
    "inconsistent": 1,
    "missing": 1,
    "mendelian_error_rate": 0.3333333333333333,
    "inconsistent_indices": [2],
    "sites": [
      {
        "index": 0,
        "site": "rs1",
        "status": "consistent",
        "paternal_allele": 0,
        "maternal_allele": 1
      },
      {
        "index": 1,
        "site": "rs2",
        "status": "consistent",
        "paternal_allele": 1,
        "maternal_allele": 1
      },
      {
        "index": 2,
        "site": "rs3",
        "status": "inconsistent"
      },
      {
        "index": 3,
        "site": "rs4",
        "status": "missing"
      }
    ],
    "rule": "a site is consistent if the child's two alleles can be formed by one allele from each parent"
  }
}
```

## 2) calculate_pedigree_likelihood

- **Purpose:** Computes the likelihood of observed biallelic genotypes on a small pedigree by exact enumeration, without genotyping error.

- **Method:** Founders follow Hardy-Weinberg proportions with alt-allele frequency p (`(1-p)^2`, `2p(1-p)`, `p^2`); a child's genotype follows Mendelian transmission, each parent passing the alt allele with probability dosage/2; sites are independent. Members that are absent from `genotypes` or have a missing call are summed out exactly by enumerating their possible genotypes (at most 8 ungenotyped members per site). The likelihood is the product over sites of the per-site probabilities; if any site has probability 0 the data are impossible and `log10_likelihood` is null. Only biallelic (0/1) genotypes are supported; no linkage between sites, no mutation, no genotyping error.

- **API:** `dnakit.haplotypes.calculate_pedigree_likelihood(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input` (required): `{"pedigree": [{"id", "father", "mother"}, ...], "genotypes": {id: [calls per site]}}`. Give both parents or neither (founders). Calls are `"0/1"`-style strings with 0/1 alleles, pairs, or dosages 0/1/2; all lists have the same length. `options` (optional): `allele_frequency` (number in (0, 1) or one number per site, default 0.5).

- **Output:** Standard envelope: `status` (`"ok"`), `module`, `operation`, `provenance`, `input_summary`, `warnings` and `result`. `result` contains `site_likelihoods`, `likelihood`, `log10_likelihood` (null if impossible), `impossible`, `individuals`, `sites` and a `model` string.

- **Example code:**

```python
import json

from dnakit.haplotypes import calculate_pedigree_likelihood

pedigree = [
    {"id": "f"},
    {"id": "m"},
    {"id": "c", "father": "f", "mother": "m"},
]
genotypes = {
    "f": ["0/0", "0/1"],
    "m": ["1/1", "0/1"],
    "c": ["0/1", "0/1"],
}
result = calculate_pedigree_likelihood(
    {"pedigree": pedigree, "genotypes": genotypes},
    options={"allele_frequency": 0.2},
)
print(json.dumps(result, indent=2))
```

- **Example output:**

```json
{
  "module": "haplotypes",
  "operation": "calculate_pedigree_likelihood",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 1,
    "text_length": 161,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "62ffa5c42689f1bfee4521017ab68ee6092ee94bf57e8516c57147ec51c9f566"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "site_likelihoods": [0.02560000000000001, 0.05120000000000002],
    "likelihood": 0.0013107200000000011,
    "log10_likelihood": -2.882490073712319,
    "impossible": false,
    "individuals": 3,
    "sites": 2,
    "model": "founders in Hardy-Weinberg proportions with the given alt-allele frequency, Mendelian transmission, no genotyping error, sites independent; ungenotyped members summed out exactly"
  }
}
```

## 3) cluster

- **Purpose:** Clusters haplotypes by single linkage on the Hamming distance (connected components of haplotypes within `max_distance` mismatches).

- **Method:** Computes the Hamming distance between every pair of haplotypes (number of differing sites; the characters `?`, `N`, `.` and `-` are skipped for a pair when `ignore_missing` is true) and joins every pair with distance <= `max_distance` with union-find; each connected component is a cluster. Single linkage means a cluster can contain haplotypes further apart than `max_distance` through a chain of close neighbours. Clusters are listed in the order of their first member. At most 5,000 haplotypes.

- **API:** `dnakit.haplotypes.cluster(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input` (required): a list of haplotype strings of equal length, or `{name: haplotype | [haplotypes]}`. `options` (optional): `max_distance` (int >= 0, default 1), `ignore_missing` (bool, default True).

- **Output:** Standard envelope: `status` (`"ok"`), `module`, `operation`, `provenance`, `input_summary`, `warnings` and `result`. `result` contains `cluster_count`, `clusters` (each with `size`, `members` as names, the distinct `haplotypes` and the `indices` of the input haplotypes), `max_distance` and a `linkage` description.

- **Example code:**

```python
import json

from dnakit.haplotypes import cluster

result = cluster(["0000", "0001", "1111", "1110"], options={"max_distance": 1})
print(json.dumps(result, indent=2))
```

- **Example output:**

```json
{
  "module": "haplotypes",
  "operation": "cluster",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 4,
    "text_length": 32,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "35c2ec48596e2606e14353d91ff34b1fa688576c2075122e1582d191f9a9069c"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "cluster_count": 2,
    "clusters": [
      {
        "size": 2,
        "members": ["hap_1", "hap_2"],
        "haplotypes": ["0000", "0001"],
        "indices": [0, 1]
      },
      {
        "size": 2,
        "members": ["hap_3", "hap_4"],
        "haplotypes": ["1110", "1111"],
        "indices": [2, 3]
      }
    ],
    "max_distance": 1,
    "linkage": "single (connected components of the graph joining haplotypes within max_distance mismatches)"
  }
}
```

## 4) compare

- **Purpose:** Compares haplotypes by pairwise Hamming distance, or compares two phased haplotype pairs by switch error rate.

- **Method:** `mode="hamming"` (default): the matrix of pairwise numbers of differing sites (`?`, `N`, `.`, `-` skipped for a pair unless `ignore_missing` is false); for exactly two haplotypes it also lists the differing positions, the distance and the number of compared sites. A sample with several haplotypes gets the names `name_1`, `name_2`, .... `mode="switch_error"`: among sites that are heterozygous in both the truth and the inferred pair and carry the same two alleles in both, the relative orientation of the two inferred haplotypes against the truth is recorded; a switch is a change of orientation between consecutive such sites and `switch_error_rate = switches / (sites - 1)` (null with fewer than two such sites). At most 2,000 haplotypes in Hamming mode.

- **API:** `dnakit.haplotypes.compare(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input` (required): for `hamming`, a list of at least two haplotype strings of equal length or `{name: haplotype | [haplotypes]}`; for `switch_error`, `{"truth": [h1, h2], "inferred": [h1, h2]}`. `options` (optional): `mode` (`"hamming"` default or `"switch_error"`), `ignore_missing` (bool, default True; Hamming mode).

- **Output:** Standard envelope: `status` (`"ok"`), `module`, `operation`, `provenance`, `input_summary`, `warnings` and `result`. `result` contains Hamming mode: `mode`, `names`, `distance_matrix`, `ignore_missing` and, for two haplotypes, `differing_positions`, `distance`, `compared_sites`. Switch-error mode: `mode`, `heterozygous_sites`, `switches`, `switch_error_rate` and a `definition` string.

- **Example code:**

```python
import json

from dnakit.haplotypes import compare

result = compare({"truth": ["0000", "1111"], "inferred": ["0011", "1100"]}, options={"mode": "switch_error"})
print(json.dumps(result, indent=2))
```

- **Example output:**

```json
{
  "module": "haplotypes",
  "operation": "compare",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 1,
    "text_length": 57,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "3268637c7959c932534b3c1dee8c926c2330c85dcc4a5353cd07a2a734bda330"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "mode": "switch_error",
    "heterozygous_sites": 4,
    "switches": 1,
    "switch_error_rate": 0.3333333333333333,
    "definition": "among heterozygous sites called identically in truth and inferred, a switch is a change of relative phase between consecutive sites; rate = switches / (sites - 1)"
  }
}
```

## 5) construct

- **Purpose:** Constructs haplotype blocks from phased 0/1 haplotypes with a complete pairwise linkage-disequilibrium rule.

- **Method:** Monomorphic sites and sites with minor-allele frequency below `min_maf` are excluded. Pairwise LD between two sites is computed from the haplotype frequencies of the given chromosomes: D = p(AB) - p(A)p(B), r^2 = D^2 / (p(A)(1-p(A))p(B)(1-p(B))) and |D'| = |D| / Dmax. Going left to right, a block is greedily extended while EVERY pair of sites in it has the chosen LD statistic >= `threshold`, and it is kept if it has at least `min_snps` sites; the next block starts after it. This is not the Gabriel et al. confidence-interval method, and D' is inflated for small samples.

- **API:** `dnakit.haplotypes.construct(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input` (required): a list of phased haplotypes (equal-length strings or lists of 0/1 alleles; at least two haplotypes, at most 2,000 sites) or `{name: haplotype(s)}`. `options` (optional): `statistic` (`"dprime"` default or `"r2"`), `threshold` (default 0.8), `min_snps` (int >= 2, default 2), `min_maf` (default 0.0), `positions` (one coordinate per site, echoed as `start_position` / `end_position` of each block).

- **Output:** Standard envelope: `status` (`"ok"`), `module`, `operation`, `provenance`, `input_summary`, `warnings` and `result`. `result` contains `blocks` (each with `first_site`, `last_site`, `sites`, `snp_count` and, if `positions` is given, `start_position` / `end_position`), `block_count`, `statistic`, `threshold`, `sites_total`, `sites_analysed`, `sites_excluded_monomorphic_or_low_maf` and a `rule` string.

- **Example code:**

```python
import json

from dnakit.haplotypes import construct

haplotypes = ["0000", "0011", "1100", "1111"]
result = construct(haplotypes, options={"positions": [100, 200, 300, 400]})
print(json.dumps(result, indent=2))
```

- **Example output:**

```json
{
  "module": "haplotypes",
  "operation": "construct",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 4,
    "text_length": 32,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "095e1dbb5d3cd7bea5f10e204f15a9aeaff2bbb73b167c7c21679c30b16c01b7"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "blocks": [
      {
        "first_site": 0,
        "last_site": 1,
        "sites": [0, 1],
        "snp_count": 2,
        "start_position": 100,
        "end_position": 200
      },
      {
        "first_site": 2,
        "last_site": 3,
        "sites": [2, 3],
        "snp_count": 2,
        "start_position": 300,
        "end_position": 400
      }
    ],
    "block_count": 2,
    "statistic": "dprime",
    "threshold": 0.8,
    "sites_total": 4,
    "sites_analysed": 4,
    "sites_excluded_monomorphic_or_low_maf": [],
    "rule": "greedy left-to-right: a block is a maximal run of >= min_snps polymorphic sites in which every pair has LD (|D'| or r^2) >= threshold; pairwise LD from the haplotype frequencies of the given chromosomes"
  }
}
```

## 6) count

- **Purpose:** Count the records of the input (non-blank lines of text, or items of a list); it does not interpret them as haplotypes.

- **Method:** Counts records in the shared local runtime (`run_candidate_operation`, helper `_local_result` in `src/dnakit/_candidate_runtime.py`): a text (or bytes, or the content of a `Path`) is split into non-blank lines and each line is one record; a list counts its items; a mapping that has one of the keys `records`, `items`, `values`, `data`, `sequences` or `reads` counts that value in the same way, any other mapping is one record. Nothing is parsed or validated: a FASTA with two sequences counts 4 (header and sequence lines), and a tab-separated line is one record whatever its content. `None` counts as one record and an empty string as zero.

- **API:** `dnakit.haplotypes.count(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input`: text, `bytes`, a `Path` of an existing file (first 1,000,000 bytes), a list, or a mapping (see Method). `options` is accepted but unused.

- **Output:** The standard envelope: `status` (`"ok"`), `module` (`"haplotypes"`), `operation` (`"count"`), `provenance` (`implementation` `"local_candidate_adapter"`, `backend`), `input_summary` (generic text summary: `record_count` = non-blank lines, `text_length`, `sequence_bases` = number of letters, `lengths`, `sha256`), `warnings` (empty) and `result`. `result` is an integer, the number of records.

- **Example code:**

```python
import json
from dnakit.haplotypes import count

result = count("H1\tACGT\nH2\tACGA\nH3\tTCGA\n")
print(json.dumps(result, indent=2))
```

- **Example output:**

```json
{
  "module": "haplotypes",
  "operation": "count",
  "provenance": {
    "implementation": "local_candidate_adapter",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 3,
    "text_length": 24,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "27781814e3c6429e6ce9b19cf327192ea6a3cb3e4364cf11bac5c2b35c79e069"
  },
  "warnings": [],
  "status": "ok",
  "result": 3
}
```

## 7) detect_compound_heterozygotes

- **Purpose:** Finds genes that carry heterozygous variants on both parental haplotypes (compound heterozygotes, variants in trans).

- **Method:** Only sites where the child is heterozygous are used. Trio mode (father and mother given): the origin of each heterozygous child allele comes from the unique Mendelian assignment of the child's alleles to (father, mother); a non-reference allele on the paternal side is a paternal variant, on the maternal side a maternal variant; sites with no unique assignment (for example both parents heterozygous) are listed as `undetermined`. Phased mode (no parents): the child's phased calls such as `0|1` decide haplotype 1 versus 2; an unphased heterozygous call raises an error. A gene is reported when it has at least `min_variants_per_haplotype` non-reference variants on each of the two haplotypes. There is no filtering by consequence, frequency or quality.

- **API:** `dnakit.haplotypes.detect_compound_heterozygotes(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input` (required): `{"genes": [gene name per site], "child": [...]}` plus either `"father"` and `"mother"` (trio mode) or no parents with PHASED child calls like `"0|1"` (phased mode). A call is `"a/b"` or `"a|b"` with integer alleles, a pair `(a, b)`, an alt-allele dosage 0/1/2, or `None` / `"./."` for missing. `options` (optional): `min_variants_per_haplotype` (int, default 1).

- **Output:** Standard envelope: `status` (`"ok"`), `module`, `operation`, `provenance`, `input_summary`, `warnings` and `result`. `result` contains `mode` (`"trio"` or `"phased_child"`), `haplotype_labels`, `genes_with_het_variants`, `compound_heterozygous_genes` (each with `gene` and the site indices `paternal`, `maternal`, `undetermined`; for phased mode `paternal` / `maternal` mean haplotype 1 / 2), `compound_heterozygous_count`, `undetermined_origin_sites` and a `rule` string.

- **Example code:**

```python
import json

from dnakit.haplotypes import detect_compound_heterozygotes

trio = {
    "genes": ["G1", "G1", "G2"],
    "child": ["0/1", "0/1", "0/1"],
    "father": ["0/1", "0/0", "0/1"],
    "mother": ["0/0", "0/1", "0/1"],
}
result = detect_compound_heterozygotes(trio)
print(json.dumps(result, indent=2))
```

- **Example output:**

```json
{
  "module": "haplotypes",
  "operation": "detect_compound_heterozygotes",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 1,
    "text_length": 127,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "0d3afcdd898ed1cdbda697dc7b7022a7b3492b561ceb8e9561c8618bd1e12284"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "mode": "trio",
    "haplotype_labels": "paternal/maternal",
    "genes_with_het_variants": 2,
    "compound_heterozygous_genes": [
      {
        "gene": "G1",
        "paternal": [0],
        "maternal": [1],
        "undetermined": []
      }
    ],
    "compound_heterozygous_count": 1,
    "undetermined_origin_sites": 1,
    "rule": "compound heterozygote = a gene with a heterozygous variant on each of the two parental haplotypes (in trans); trio mode assigns origin when only one parent carries the alt allele, phased mode uses the phase of the child's calls"
  }
}
```

## 8) detect_denovo_variants

- **Purpose:** Detects candidate de novo variants in a trio: child alleles that are absent from both parents' genotypes.

- **Method:** A site is a de novo candidate when the child carries an allele that is contained in neither the father's nor the mother's genotype (for example parents 0/0 and 0/0, child 0/1). Other Mendelian violations (for example a parental allele that is missing from the child's genotype) are listed separately as `other_mendelian_inconsistent_indices`. A site with any missing call is skipped. Read depth, genotype quality and allele balance are not available, so no such filtering is applied; this is a genotype-level screen, not a de novo caller.

- **API:** `dnakit.haplotypes.detect_denovo_variants(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input` (required): `{"child": [...], "father": [...], "mother": [...]}` with one genotype call per site (lists of equal length) and an optional `"sites"` list of labels. A call is `"a/b"` or `"a|b"` with integer alleles, a pair `(a, b)`, an alt-allele dosage 0/1/2, or `None` / `"./."` for missing. `options`: none.

- **Output:** Standard envelope: `status` (`"ok"`), `module`, `operation`, `provenance`, `input_summary`, `warnings` and `result`. `result` contains `evaluated_sites`, `denovo_count`, `denovo` (each with `index`, optional `site` label, `child_genotype` and `novel_alleles`), `other_mendelian_inconsistent_indices` and a `rule` string.

- **Example code:**

```python
import json

from dnakit.haplotypes import detect_denovo_variants

trio = {
    "child": ["0/1", "0/1", "1/1", "0/0"],
    "father": ["0/0", "0/1", "0/0", "0/0"],
    "mother": ["0/0", "0/0", "0/1", "0/0"],
}
result = detect_denovo_variants(trio)
print(json.dumps(result, indent=2))
```

- **Example output:**

```json
{
  "module": "haplotypes",
  "operation": "detect_denovo_variants",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 1,
    "text_length": 119,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "60c25594b0884721304b995108bca1551648ed3b3b4291db0bd284e5c721c7e7"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "evaluated_sites": 4,
    "denovo_count": 1,
    "denovo": [
      {
        "index": 0,
        "child_genotype": "0/1",
        "novel_alleles": [1]
      }
    ],
    "other_mendelian_inconsistent_indices": [2],
    "rule": "de novo = the child carries an allele that neither parent's genotype contains; other Mendelian violations (e.g. missing parental allele transmitted) are listed separately"
  }
}
```

## 9) frequency

- **Purpose:** Counts haplotypes and computes their frequencies and Nei's haplotype diversity.

- **Method:** Counts the distinct haplotype strings. frequency = count / total number of chromosomes; haplotype diversity (Nei) = n/(n-1) * (1 - sum p_i^2), null when n < 2. Rows are sorted by decreasing count, then by haplotype string. A mapping whose values are all integers is read as `{haplotype: count}`.

- **API:** `dnakit.haplotypes.frequency(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input` (required): a list of haplotype strings, `{sample: [h1, h2]}`, `{name: haplotype}` or `{haplotype: count}`. `options`: none.

- **Output:** Standard envelope: `status` (`"ok"`), `module`, `operation`, `provenance`, `input_summary`, `warnings` and `result`. `result` contains `total` (chromosomes), `distinct`, `haplotypes` (rows with `haplotype`, `count`, `frequency`), `haplotype_diversity` and a `definition` string.

- **Example code:**

```python
import json

from dnakit.haplotypes import frequency

result = frequency({"s1": ["01", "10"], "s2": ["01", "01"], "s3": ["11", "01"]})
print(json.dumps(result, indent=2))
```

- **Example output:**

```json
{
  "module": "haplotypes",
  "operation": "frequency",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 1,
    "text_length": 60,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "d5e0fa610438f0870bd6a844d434431de3eb2c92f255ac64c3012dffc59bd5d1"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "total": 6,
    "distinct": 3,
    "haplotypes": [
      {
        "haplotype": "01",
        "count": 4,
        "frequency": 0.6666666666666666
      },
      {
        "haplotype": "10",
        "count": 1,
        "frequency": 0.16666666666666666
      },
      {
        "haplotype": "11",
        "count": 1,
        "frequency": 0.16666666666666666
      }
    ],
    "haplotype_diversity": 0.6,
    "definition": "frequency = count / total chromosomes; haplotype diversity (Nei) = n/(n-1) * (1 - sum p_i^2)"
  }
}
```

## 10) impute

- **Purpose:** Genotype imputation is not implemented: it needs a phased reference haplotype panel and an HMM imputation engine.

- **Method:** Imputation infers missing genotypes from a reference panel with a statistical (HMM) model. Neither the panel nor such a model is part of the input, and no heuristic stand-in is applied, so the function never fabricates genotypes. It returns a structured `not_implemented` envelope for any input. Use Beagle, IMPUTE5 or Minimac4 outside DNAKit.

- **API:** `dnakit.haplotypes.impute(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input` (required): anything (it is only summarised in `input_summary`). `options` (optional): ignored.

- **Output:** Always returns `status: "not_implemented"` (with `module`, `operation`, `provenance`, `input_summary`, `warnings`), `required_dependency` (what is missing: a phased reference haplotype panel and an imputation backend such as Beagle, IMPUTE5 or Minimac4) and `message` (why it cannot be computed from the input). No `result` is produced.

- **Example code:**

```python
import json

from dnakit.haplotypes import impute

result = impute({"sample1": ["0/1", "./.", "1/1"]})
print(json.dumps(result, indent=2))
```

- **Example output:**

```json
{
  "module": "haplotypes",
  "operation": "impute",
  "provenance": {
    "implementation": "not_implemented",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 1,
    "text_length": 34,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "89df9963733624e389004dc0717e7ddaf1bb4794e04960464957ab1395e8e1c2"
  },
  "warnings": [],
  "status": "not_implemented",
  "required_dependency": "phased reference haplotype panel and an imputation backend (Beagle, IMPUTE5, Minimac4)",
  "message": "Imputation infers missing genotypes from a reference panel with a statistical model; neither is part of the input."
}
```

## 11) infer_parentage

- **Purpose:** Screens candidate parents for a child by Mendelian exclusion (sites where parent and child share no allele).

- **Method:** For every candidate, compares the genotype at each site where both child and candidate are called and counts the sites at which they share no allele (for example opposite homozygotes 0/0 vs 1/1), since a true parent must transmit one allele to the child. A candidate is excluded when this count exceeds `max_mismatches`. Among the non-excluded candidates the ones with the fewest mismatches are the best candidates. This is exclusion only, not a paternity or kinship likelihood, and it uses no allele frequencies or genotyping-error model. A candidate with no jointly called site cannot be assessed: it is neither excluded nor listed in `not_excluded` or `best_candidates`, but is reported in `not_assessed` and named in a warning.

- **API:** `dnakit.haplotypes.infer_parentage(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input` (required): `{"child": [calls], "candidates": {name: [calls]}}`, one call per site and all lists of equal length. A call is `"a/b"` or `"a|b"` with integer alleles, a pair `(a, b)`, an alt-allele dosage 0/1/2, or `None` / `"./."` for missing. `options` (optional): `max_mismatches` (int >= 0, default 0).

- **Output:** Standard envelope: `status` (`"ok"`), `module`, `operation`, `provenance`, `input_summary`, `warnings` and `result`. `result` contains `candidates` (per candidate `compared_sites`, `mismatches`, `mismatch_rate`, `excluded`), `not_excluded`, `not_assessed` (candidate ids with no jointly called site, sorted), `best_candidates`, `max_mismatches` and a `rule` string; the envelope `warnings` lists the candidates that could not be assessed.

- **Example code:**

```python
import json

from dnakit.haplotypes import infer_parentage

data = {
    "child": ["0/0", "1/1", "0/1"],
    "candidates": {
        "A": ["0/1", "1/1", "0/0"],
        "B": ["1/1", "0/0", "0/1"],
    },
}
result = infer_parentage(data)
print(json.dumps(result, indent=2))
```

- **Example output:**

```json
{
  "module": "haplotypes",
  "operation": "infer_parentage",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 1,
    "text_length": 104,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "501c4881df823b8b5ba158820fe055215be69bcb17b1bcb9ba0d68f4d3085860"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "candidates": {
      "A": {
        "compared_sites": 3,
        "mismatches": 0,
        "mismatch_rate": 0.0,
        "excluded": false
      },
      "B": {
        "compared_sites": 3,
        "mismatches": 2,
        "mismatch_rate": 0.6666666666666666,
        "excluded": true
      }
    },
    "not_excluded": ["A"],
    "not_assessed": [],
    "best_candidates": ["A"],
    "max_mismatches": 0,
    "rule": "exclusion method: a candidate parent is excluded when the number of sites where parent and child share no allele (e.g. opposite homozygotes) exceeds max_mismatches; a candidate with no jointly called site cannot be assessed and is listed in not_assessed (never in not_excluded or best_candidates)"
  }
}
```

## 12) phase

- **Purpose:** Phases a child's genotypes with Mendelian trio logic (paternal|maternal order).

- **Method:** A heterozygous child site is phased when exactly one assignment of its two alleles to (father, mother) is compatible with the parents' genotypes; the call is then written `paternal|maternal`. If both parents are heterozygous (or otherwise more than one assignment fits) the site stays unphased (`a/b`, status `ambiguous`). Mendelian violations are flagged (`mendelian_inconsistent`); a homozygous child site is written `a|a`; a missing child call stays `None`; if a parent is missing at a heterozygous child site the status is `unresolved_parent_missing`. Statistical (population) phasing and read-based phasing are not performed, so phase is not propagated between sites.

- **API:** `dnakit.haplotypes.phase(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input` (required): `{"child": [...], "father": [...], "mother": [...]}` with one genotype call per site (lists of equal length) and an optional `"sites"` list of labels. A call is `"a/b"` or `"a|b"` with integer alleles, a pair `(a, b)`, an alt-allele dosage 0/1/2, or `None` / `"./."` for missing. `options`: none.

- **Output:** Standard envelope: `status` (`"ok"`), `module`, `operation`, `provenance`, `input_summary`, `warnings` and `result`. `result` contains `phased_child` (one call per site: `"a|b"` phased, `"a/b"` unphased, or null), `status` (per site: `phased`, `ambiguous`, `homozygous`, `missing`, `mendelian_inconsistent` or `unresolved_parent_missing`), `phased_het_sites`, `ambiguous_het_sites`, `inconsistent_sites`, `fraction_het_sites_phased`, `convention` (`"paternal|maternal"`) and a `rule` string.

- **Example code:**

```python
import json

from dnakit.haplotypes import phase

trio = {
    "child": ["0/1", "0/1", "0/1", "1/1", "./.", "0/1"],
    "father": ["0/0", "0/1", "0/0", "1/1", "0/0", "0/1"],
    "mother": ["0/1", "0/1", "1/1", "0/1", "0/0", "0/0"],
}
result = phase(trio)
print(json.dumps(result, indent=2))
```

- **Example output:**

```json
{
  "module": "haplotypes",
  "operation": "phase",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 1,
    "text_length": 161,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "81ac1fa8bb2f1d1a55174a56e17a30af2f26a94ebd9aa8ead62f94a35869ad8f"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "phased_child": ["0|1", "0/1", "0|1", "1|1", null, "1|0"],
    "status": ["phased", "ambiguous", "phased", "homozygous", "missing", "phased"],
    "phased_het_sites": 3,
    "ambiguous_het_sites": 1,
    "inconsistent_sites": 0,
    "fraction_het_sites_phased": 0.75,
    "convention": "paternal|maternal",
    "rule": "Mendelian trio phasing: a heterozygous child site is phased when exactly one assignment of its two alleles to (father, mother) is compatible with the parents' genotypes"
  }
}
```

## 13) reconstruct

- **Purpose:** Estimates haplotype frequencies and each individual's haplotype pair from unphased biallelic genotypes by the EM algorithm.

- **Method:** Excoffier-Slatkin (1995) EM under Hardy-Weinberg equilibrium. For each individual all haplotype pairs compatible with its genotype are enumerated (at most 12 heterozygous sites per individual, at most 4,096 candidate haplotypes in total). Starting from uniform haplotype frequencies, the E step weights each pair by `freq(h1) * freq(h2)` (times 2 for heterozygous pairs) and the M step re-estimates the frequencies from the expected chromosome counts, until the largest change in a frequency is below `tolerance` or `max_iterations` is reached. Each individual is assigned its maximum-posterior haplotype pair. The reported `log_likelihood` is the log-likelihood of the data at the final haplotype frequencies (computed after the last M step), so it matches the returned `frequencies`. Missing genotypes and multi-allelic sites are not supported, and there is no recombination or linkage model.

- **API:** `dnakit.haplotypes.reconstruct(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input` (required): `{sample: [calls per site]}` with calls as alt dosage 0/1/2 or `"0/1"`-style strings, no missing calls, the same number of sites for every sample. `options` (optional): `max_iterations` (int, default 1000), `tolerance` (maximum frequency change to stop, default 1e-8).

- **Output:** Standard envelope: `status` (`"ok"`), `module`, `operation`, `provenance`, `input_summary`, `warnings` and `result`. `result` contains `frequencies` (haplotypes with estimated frequency above 1e-12, sorted by decreasing frequency), `individuals` (per sample: the most probable `haplotypes` pair, its `posterior` and the number of `alternatives`), `iterations`, `converged`, `log_likelihood` and a `model` string.

- **Example code:**

```python
import json

from dnakit.haplotypes import reconstruct

genotypes = {
    "a": [0, 0],
    "b": [0, 0],
    "c": [2, 2],
    "d": [2, 2],
    "e": [1, 1],
}
result = reconstruct(genotypes)
print(json.dumps(result, indent=2))
```

- **Example output:**

```json
{
  "module": "haplotypes",
  "operation": "reconstruct",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 1,
    "text_length": 65,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "0f79d304afa9d889f84bd34456e47910635e840ace689d24058c1d3a3d2b5d85"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "frequencies": [
      {
        "haplotype": "00",
        "frequency": 0.5
      },
      {
        "haplotype": "11",
        "frequency": 0.5
      }
    ],
    "individuals": {
      "a": {
        "haplotypes": ["00", "00"],
        "posterior": 1.0,
        "alternatives": 1
      },
      "b": {
        "haplotypes": ["00", "00"],
        "posterior": 1.0,
        "alternatives": 1
      },
      "c": {
        "haplotypes": ["11", "11"],
        "posterior": 1.0,
        "alternatives": 1
      },
      "d": {
        "haplotypes": ["11", "11"],
        "posterior": 1.0,
        "alternatives": 1
      },
      "e": {
        "haplotypes": ["00", "11"],
        "posterior": 1.0,
        "alternatives": 2
      }
    },
    "iterations": 5,
    "converged": true,
    "log_likelihood": -6.238324625039508,
    "model": "Excoffier-Slatkin (1995) EM for haplotype frequencies from unphased biallelic genotypes under Hardy-Weinberg, uniform initial frequencies; individual haplotypes = maximum posterior pair"
  }
}
```

## 14) select

- **Purpose:** Selects haplotypes by minimum frequency, minimum count and/or top-N from a haplotype sample.

- **Method:** Builds the frequency table as in `frequency` (frequency = count / total chromosomes), keeps the haplotypes with frequency >= `min_frequency` and count >= `min_count`, then truncates to the first `top_n` rows of the table sorted by decreasing frequency. `covered_frequency` is the summed frequency of the retained haplotypes. No LD or functional criterion is used.

- **API:** `dnakit.haplotypes.select(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input` (required): as for `frequency` (a list of haplotype strings, `{sample: [h1, h2]}` or `{haplotype: count}`). `options` (optional): `min_frequency` (default 0.05), `min_count` (int, default 1), `top_n` (int, default none).

- **Output:** Standard envelope: `status` (`"ok"`), `module`, `operation`, `provenance`, `input_summary`, `warnings` and `result`. `result` contains `selected` (rows with `haplotype`, `count`, `frequency`, sorted by decreasing frequency), `selected_count`, `total_distinct`, `covered_frequency`, `min_frequency`, `min_count` and `top_n`.

- **Example code:**

```python
import json

from dnakit.haplotypes import select

sample = ["01"] * 5 + ["10"] * 4 + ["11"]
result = select(sample, options={"min_frequency": 0.15})
print(json.dumps(result, indent=2))
```

- **Example output:**

```json
{
  "module": "haplotypes",
  "operation": "select",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 10,
    "text_length": 60,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "e80a1b0d3639709d7057db6f0431259c8919ad8b80b9dce7ee4d93f1f54f67fb"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "selected": [
      {
        "haplotype": "01",
        "count": 5,
        "frequency": 0.5
      },
      {
        "haplotype": "10",
        "count": 4,
        "frequency": 0.4
      }
    ],
    "selected_count": 2,
    "total_distinct": 3,
    "covered_frequency": 0.9,
    "min_frequency": 0.15,
    "min_count": 1,
    "top_n": null
  }
}
```

## 15) simulate

- **Purpose:** Simulates diploid individuals by drawing two haplotypes from given haplotype frequencies (random mating).

- **Method:** Each of the two chromosomes of each individual is drawn independently, with replacement, from the haplotypes weighted by the given frequencies or counts (normalised to sum to 1), using Python's `random.Random(seed)`; the same `seed` reproduces the result. No mutation, recombination or demography is simulated; this is not a coalescent simulator. Haplotype strings are only labels and are not interpreted.

- **API:** `dnakit.haplotypes.simulate(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input` (required): `{haplotype: frequency or count}` (non-negative numbers, not all zero). `options` (optional): `n_individuals` (int 1-100,000, default 10), `seed` (int, default 0).

- **Output:** Standard envelope: `status` (`"ok"`), `module`, `operation`, `provenance`, `input_summary`, `warnings` and `result`. `result` contains `individuals` (each with `id` `sim_<k>` and its two `haplotypes`), `observed_haplotype_counts`, `input_frequencies` (normalised), `n_individuals`, `seed` and a `model` string.

- **Example code:**

```python
import json

from dnakit.haplotypes import simulate

result = simulate({"01": 0.3, "10": 0.7}, options={"n_individuals": 3, "seed": 7})
print(json.dumps(result, indent=2))
```

- **Example output:**

```json
{
  "module": "haplotypes",
  "operation": "simulate",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 1,
    "text_length": 22,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "f1a0a6b72feeb29090802690be6ea8477c7f021b831c1a860459360b5dd23bfd"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "individuals": [
      {
        "id": "sim_1",
        "haplotypes": ["10", "01"]
      },
      {
        "id": "sim_2",
        "haplotypes": ["10", "01"]
      },
      {
        "id": "sim_3",
        "haplotypes": ["10", "10"]
      }
    ],
    "observed_haplotype_counts": {
      "01": 2,
      "10": 4
    },
    "input_frequencies": {
      "01": 0.3,
      "10": 0.7
    },
    "n_individuals": 3,
    "seed": 7,
    "model": "random mating: each of the 2 chromosomes of an individual is drawn independently from the haplotype frequencies (python random.Random(seed)); no mutation or recombination"
  }
}
```

## 16) tag

- **Purpose:** Selects tag SNPs from phased 0/1 haplotypes by greedy r^2 binning (Carlson et al. 2004).

- **Method:** Monomorphic sites and sites with minor-allele frequency below `min_maf` are excluded. Pairwise r^2 is computed from the haplotype frequencies of the given chromosomes. Repeatedly, the site whose r^2 with the largest number of still-untagged sites (itself included) is >= `r2_threshold` becomes a tag (ties: lowest index) and those sites form its bin and become tagged, until every site is tagged. Every site in a bin has r^2 >= threshold with the tag, but not necessarily with the other sites of the bin.

- **API:** `dnakit.haplotypes.tag(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input` (required): a list of phased haplotypes (equal-length strings or lists of 0/1 alleles; at least two haplotypes, at most 2,000 sites) or `{name: haplotype(s)}`. `options` (optional): `r2_threshold` (default 0.8), `min_maf` (default 0.0).

- **Output:** Standard envelope: `status` (`"ok"`), `module`, `operation`, `provenance`, `input_summary`, `warnings` and `result`. `result` contains `tag_sites` (site indices), `bins` (tag site index -> sites it tags), `tag_count`, `sites_analysed`, `r2_threshold`, `excluded_monomorphic_or_low_maf` and a `rule` string.

- **Example code:**

```python
import json

from dnakit.haplotypes import tag

result = tag(["0000", "0011", "1100", "1111"], options={"r2_threshold": 0.8})
print(json.dumps(result, indent=2))
```

- **Example output:**

```json
{
  "module": "haplotypes",
  "operation": "tag",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 4,
    "text_length": 32,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "095e1dbb5d3cd7bea5f10e204f15a9aeaff2bbb73b167c7c21679c30b16c01b7"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "tag_sites": [0, 2],
    "bins": {
      "0": [0, 1],
      "2": [2, 3]
    },
    "tag_count": 2,
    "sites_analysed": 4,
    "r2_threshold": 0.8,
    "excluded_monomorphic_or_low_maf": [],
    "rule": "greedy tagging (Carlson et al. 2004): repeatedly choose the site with r^2 >= threshold to the most still-untagged sites (ties: lowest index), assign those sites to its bin"
  }
}
```

## 17) visualize

- **Purpose:** Draws haplotypes as an SVG allele grid (rows = haplotypes, columns = sites).

- **Method:** One labelled row per haplotype and one square per site; allele `1` is dark red, allele `0` light blue-grey and any other character (for example `N`) grey. Each square carries a tooltip `<title>` with the haplotype, site index and allele. No reordering or clustering is done. At most 200,000 cells; pure Python, no plotting library. If `output` is given the SVG is also written to that file.

- **API:** `dnakit.haplotypes.visualize(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input` (required): a list of haplotype strings of equal length, or `{name: haplotype | [haplotypes]}`. `options` (optional): `output` (file path for the SVG; default none).

- **Output:** Standard envelope: `status` (`"ok"`), `module`, `operation`, `provenance`, `input_summary`, `warnings` and `result`. `result` contains `svg` (the SVG text), `haplotype_count`, `site_count` and `output_path` (the file written, or null).

- **Example code:**

```python
import json

from dnakit.haplotypes import visualize

result = visualize({"h1": "01", "h2": "10", "h3": "1N"})
print(json.dumps(result, indent=2))
```

- **Example output:**

```json
{
  "module": "haplotypes",
  "operation": "visualize",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 1,
    "text_length": 36,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "f161d5549dd126e6deff5d0160767f75f7ee3854cbc850ab2e9865de4fa0089c"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "svg": "<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"64\" height=\"50\" viewBox=\"0 0 64 50\">\n<text x=\"4\" y=\"15\" font-size=\"11\">h1</text>\n<rect x=\"28\" y=\"4\" width=\"13\" height=\"13\" fill=\"#dfe7ef\"><title>h1 site 0: 0</title></rect>\n<rect x=\"42\" y=\"4\" w...",
    "haplotype_count": 3,
    "site_count": 2,
    "output_path": null
  }
}
```

(Long arrays and text in the example output are shortened with `...`.)
