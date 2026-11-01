# Assembly, variant and GWAS command-line facades: GEMMA, QUAST, BUSCO and DELLY

This group of functions runs four caller-installed programs through fixed argument whitelists, without a shell: [GEMMA](https://github.com/genetics-statistics/GEMMA) (variance components from summary statistics, eigen decomposition, G x E, multivariate LMM, BSLMM, prediction), [QUAST](https://github.com/ablab/quast) (reference-based assembly evaluation, gene and rRNA prediction, read mapping, MetaQUAST), [BUSCO](https://busco.ezlab.org/) (completeness with a chosen gene predictor, assembly statistics, scaffold composition) and [DELLY](https://github.com/dellytools/delly) (genotyping, assembly-based SVs, copy numbers, merge, filter, breakend conversion, read-depth plot). Missing programs raise an error with the install command; output files and folders must not exist yet.

Tested versions: GEMMA 0.98.5, QUAST 5.2.0, DELLY 2.7.0. Notes: GEMMA 0.98.5 aborts with a heap error after `-gs` has written its files (the facade accepts this and reports it in `warnings`); its `-gw` (LD-score weights) is not implemented upstream and its BIMBAM `-gxe` reads the wrong file, so the G x E function takes PLINK files. BUSCO and R were not available in the test environment: the BUSCO facade is tested with a stand-in program that mimics its command line and output files, and the read-depth plot (`R/rd.R`) is documented but was not run. The example call for BUSCO and the plot are therefore not executed.

- **Example code:**

```python
import random
from pathlib import Path

import pysam

N, P = 80, 60


def revcomp(s: str) -> str:
    return s[::-1].translate(str.maketrans("ACGT", "TGCA"))


def build(folder: Path) -> dict[str, Path]:
    """Reference, short-read BAM, long-read BAM and assembly BAM that all carry one 800 bp deletion at 30000."""
    random.seed(11)
    length, start, size = 60000, 30000, 800
    ref = "".join(random.choice("ACGT") for _ in range(length))
    donor = ref[:start] + ref[start + size :]
    (folder / "ref.fa").write_text(">chr1\n" + "\n".join(ref[i : i + 70] for i in range(0, length, 70)) + "\n")
    pysam.faidx(str(folder / "ref.fa"))
    header = {"HD": {"VN": "1.6", "SO": "coordinate"}, "SQ": [{"SN": "chr1", "LN": length}], "RG": [{"ID": "s1", "SM": "sample1"}]}

    def seg(name, seq, pos, cigar, flag, tags=()):
        a = pysam.AlignedSegment()
        a.query_name, a.query_sequence, a.flag = name, seq, flag
        a.reference_id, a.reference_start, a.mapping_quality, a.cigarstring = 0, pos, 60, cigar
        a.query_qualities = pysam.qualitystring_to_array("I" * len(seq))
        a.set_tag("RG", "s1")
        for tag, value in tags:
            a.set_tag(tag, value)
        return a

    def to_ref(pos):
        return pos if pos < start else pos + size

    def write(name, records):
        records.sort(key=lambda r: (r.reference_start, r.query_name))
        with pysam.AlignmentFile(str(folder / name), "wb", header=header) as out:
            for r in records:
                out.write(r)
        pysam.index(str(folder / name))

    short = []
    for i in range(int(len(donor) * 25 / 200)):
        insert = max(200, int(random.gauss(350, 30)))
        s = random.randrange(0, len(donor) - insert)
        frag = donor[s : s + insert]
        pair = []
        for seq, dpos, flag in ((frag[:100], s, 0x1 | 0x2 | 0x20 | 0x40), (revcomp(frag)[:100], s + insert - 100, 0x1 | 0x2 | 0x10 | 0x80)):
            if dpos < start < dpos + 100:
                k = start - dpos
                pair.append(seg(f"p{i}", seq, to_ref(dpos), f"{k}M{100 - k}S", flag, [("SA", f"chr1,{start + size + 1},+,{k}S{100 - k}M,60,0;")]))
            else:
                pair.append(seg(f"p{i}", seq, to_ref(dpos), "100M", flag))
        tlen = pair[1].reference_start + 100 - pair[0].reference_start
        for x, y, t in ((pair[0], pair[1], tlen), (pair[1], pair[0], -tlen)):
            x.next_reference_id, x.next_reference_start, x.template_length = 0, y.reference_start, t
            if tlen > 700:
                x.flag &= ~0x2
        short += pair
        for x in pair:
            if x.has_tag("SA"):
                k = x.cigartuples[0][1]
                sup = seg(x.query_name, x.query_sequence, start + size, f"{k}S{100 - k}M", x.flag | 0x800, [("SA", f"chr1,{x.reference_start + 1},+,{k}M{100 - k}S,60,0;")])
                sup.next_reference_id, sup.next_reference_start, sup.template_length = 0, x.next_reference_start, x.template_length
                short.append(sup)
    write("short.bam", short)
    long_reads = []
    for i in range(20):
        s = start - random.randrange(1500, 3500)
        seq = donor[s : s + 6000]
        k = start - s
        long_reads.append(seg(f"lr{i}", seq, s, f"{k}M{len(seq) - k}S", 0, [("SA", f"chr1,{start + size + 1},+,{k}S{len(seq) - k}M,60,0;")]))
        long_reads.append(seg(f"lr{i}", seq, start + size, f"{k}S{len(seq) - k}M", 0x800, [("SA", f"chr1,{s + 1},+,{k}M{len(seq) - k}S,60,0;")]))
    write("long.bam", long_reads)
    contig = donor[start - 10000 : start + 10000]
    write("asm.bam", [
        seg("contig1", contig, start - 10000, "10000M10000S", 0, [("SA", f"chr1,{start + size + 1},+,10000S10000M,60,0;")]),
        seg("contig1", contig, start + size, "10000S10000M", 0x800, [("SA", f"chr1,{start - 9999},+,10000M10000S,60,0;")]),
    ])
    return {name: folder / name for name in ("ref.fa", "short.bam", "long.bam", "asm.bam")}


def toy(folder: Path) -> dict[str, Path]:
    random.seed(7)
    freq = [random.uniform(0.2, 0.5) for _ in range(P)]
    geno = [[sum(random.random() < f for _ in range(2)) for _ in range(N)] for f in freq]
    (folder / "geno.txt").write_text("".join(f"rs{j + 1}, A, G, " + ", ".join(map(str, row)) + "\n" for j, row in enumerate(geno)))
    rows = []
    for i in range(N):
        y1 = 0.8 * geno[3][i] - 0.6 * geno[20][i] + random.gauss(0, 1)
        rows.append((y1, 0.5 * y1 + random.gauss(0, 1)))
    (folder / "pheno.txt").write_text("".join(f"{'NA' if i >= 70 else format(a, '.4f')}\t{b:.4f}\n" for i, (a, b) in enumerate(rows)))
    (folder / "env.txt").write_text("".join(f"{random.gauss(0, 1):.4f}\n" for _ in range(N)))
    code = {2: 0b00, 1: 0b10, 0: 0b11}
    with (folder / "toy.bed").open("wb") as bed, (folder / "toy.bim").open("w") as bim:
        bed.write(bytes([0x6C, 0x1B, 0x01]))
        for j, row in enumerate(geno):
            out = bytearray()
            for k in range(0, N, 4):
                out.append(sum(code[row[k + m]] << (2 * m) for m in range(4) if k + m < N))
            bed.write(out)
            bim.write(f"1\trs{j + 1}\t0\t{1000 * (j + 1)}\tA\tG\n")
    (folder / "toy.fam").write_text("".join(f"f{i} i{i} 0 0 1 -9\n" for i in range(N)))
    (folder / "beta.txt").write_text("SNP N Z INC_ALLELE DEC_ALLELE\n" + "".join(f"rs{j + 1} 70 {((j * 7) % 11 - 5) / 5.0} A G\n" for j in range(P)))
    return {name: folder / name for name in ("geno.txt", "pheno.txt", "env.txt", "beta.txt")}


def genome_files(folder: Path) -> tuple[Path, Path]:
    random.seed(21)
    genome = "".join(random.choice("ACGT") for _ in range(20000))
    reference = folder / "ref.fa"
    reference.write_text(f">ref\n{genome}\n")
    contigs = folder / "contigs.fa"
    contigs.write_text(f">c1\n{genome[:9000]}\n>c2\n{genome[9500:]}\n")
    return reference, contigs


for name in ("sv", "gwas", "asm"):
    Path(name).mkdir(exist_ok=True)
build(Path("sv"))
toy(Path("gwas"))
genome_files(Path("asm"))
print(sorted(p.name for p in Path("sv").glob("*.bam")), len(list(Path("gwas").iterdir())), sorted(p.name for p in Path("asm").iterdir()))
```

- **Example output:**

```text
['asm.bam', 'long.bam', 'short.bam'] 7 ['contigs.fa', 'ref.fa']
```

## 1) `dnakit.gwas.gemma_p15_cli_functions.gemma_snp_covariance_matrix()`

- **Purpose:** Compute the SNP variance-covariance (S) matrix for summary-statistic variance components.
- **Method:** Runs `gemma -gs`. The S matrix of the SNPs of each category is estimated from the individuals with a phenotype. GEMMA 0.98.5 aborts after writing the files; this is accepted and reported in `warnings`.
- **API:** `dnakit.gwas.gemma_p15_cli_functions.gemma_snp_covariance_matrix(output_directory, *, plink_prefix=None, genotype_file=None, phenotype_file=None, category_file=None, sample=None, output_prefix='result', miss=None, maf=None, hwe=None, r2=None, not_snp=False, gemma_path=None)`.
- **Input:** BIMBAM or PLINK genotypes, a phenotype file, optional `category_file`, `sample`.
- **Output:** `s_matrix`, `sizes`, `warnings`, files in the new folder.
- **Example code:**

```python
import subprocess
from dnakit.gwas.gemma_p15_cli_functions import *

g, p = 'gwas/geno.txt', 'gwas/pheno.txt'
import tempfile
from pathlib import Path

out = Path(tempfile.mkdtemp())
result = gemma_snp_covariance_matrix(str(out / "s"), genotype_file=g, phenotype_file=p)
print(result["sizes"], result["status"], len(result["warnings"]))
print(round(result["s_matrix"][0][0], 5))
```

- **Example output:**

```text
[60.0, 70.0] completed 1
0.01604
```

## 2) `dnakit.gwas.gemma_p15_cli_functions.gemma_summary_statistic_q_vector()`

- **Purpose:** Compute the q vector of a GWAS summary-statistic file.
- **Method:** Runs `gemma -gq`. GEMMA 0.98.5 also needs a phenotype file to start; its values are not used.
- **API:** `dnakit.gwas.gemma_p15_cli_functions.gemma_summary_statistic_q_vector(output_directory, beta_file, phenotype_file, *, category_file=None, output_prefix='result', gemma_path=None)`.
- **Input:** A beta/z file (SNP, N, Z, alleles), a phenotype file, optional category file.
- **Output:** `q`, `q_errors`, `sizes`.
- **Example code:**

```python
import subprocess
from dnakit.gwas.gemma_p15_cli_functions import *

g, p = 'gwas/geno.txt', 'gwas/pheno.txt'
import tempfile
from pathlib import Path

out = Path(tempfile.mkdtemp())
result = gemma_summary_statistic_q_vector(str(out / "q"), "gwas/beta.txt", p)
print(result["sizes"], [round(x, 5) for x in result["q"]])
```

- **Example output:**

```text
[55.0, 70.0] [-0.00796]
```

## 3) `dnakit.gwas.gemma_p15_cli_functions.gemma_eigen_decomposition()`

- **Purpose:** Eigen-decompose a kinship (relatedness) matrix.
- **Method:** Runs `gemma -eigen`: takes the rows of the analysed individuals, centres the matrix and decomposes it.
- **API:** `dnakit.gwas.gemma_p15_cli_functions.gemma_eigen_decomposition(output_directory, kinship_file, *, plink_prefix=None, genotype_file=None, phenotype_file=None, covariates_file=None, output_prefix='result', gemma_path=None)`.
- **Input:** A kinship file from `gemma -gk`, the genotype and phenotype files, a new output folder.
- **Output:** `eigen_values`, `individual_count`; `<prefix>.eigenU.txt` holds the eigenvectors.
- **Example code:**

```python
import subprocess
from dnakit.gwas.gemma_p15_cli_functions import *

g, p = 'gwas/geno.txt', 'gwas/pheno.txt'
import tempfile
from pathlib import Path

out = Path(tempfile.mkdtemp())
subprocess.run(["gemma", "-g", g, "-p", p, "-gk", "1", "-outdir", str(out), "-o", "kin"], check=True, capture_output=True)
result = gemma_eigen_decomposition(str(out / "eig"), str(out / "kin.cXX.txt"), genotype_file=g, phenotype_file=p)
print(result["individual_count"], round(result["eigen_values"][-1], 4))
```

- **Example output:**

```text
70 1.7967
```

## 4) `dnakit.gwas.gemma_p15_cli_functions.gemma_summary_variance_components()`

- **Purpose:** Estimate variance components and heritability from GWAS summary statistics and a reference LD panel.
- **Method:** Runs `gemma -vc` with the MQS-HEW (`-vc 1`) or MQS-LDW (`-vc 2`, needs an LD-score file) algorithm.
- **API:** `dnakit.gwas.gemma_p15_cli_functions.gemma_summary_variance_components(output_directory, beta_file, *, plink_prefix=None, genotype_file=None, phenotype_file=None, category_file=None, ld_score_file=None, method='mqs_hew', covariates_file=None, sample=None, output_prefix='result', miss=None, maf=None, hwe=None, r2=None, not_snp=False, gemma_path=None)`.
- **Input:** A beta/z file, genotypes (PLINK or BIMBAM), a phenotype file, optional categories and `method`.
- **Output:** `log_values` with `pve_estimates`, `se_pve`, `sigma2_estimates`, `enrichment`.
- **Example code:**

```python
import subprocess
from dnakit.gwas.gemma_p15_cli_functions import *

g, p = 'gwas/geno.txt', 'gwas/pheno.txt'
import tempfile
from pathlib import Path

out = Path(tempfile.mkdtemp())
result = gemma_summary_variance_components(str(out / "vc"), "gwas/beta.txt", genotype_file=g, phenotype_file=p)
values = result["log_values"]
print(values["number_of_analyzed_snps_var"], round(values["pve_estimates"], 4), round(values["se_pve"], 4))
```

- **Example output:**

```text
55.0 -0.4961 0.0203
```

## 5) `dnakit.gwas.gemma_p15_cli_functions.gemma_summary_variance_component_errors()`

- **Purpose:** Standard errors of summary-statistic variance components with the asymptotic method and a reference panel.
- **Method:** Runs `gemma -ci` with `-ref` (prefix of an earlier fit) and `-pve`.
- **API:** `dnakit.gwas.gemma_p15_cli_functions.gemma_summary_variance_component_errors(output_directory, beta_file, reference_prefix, pve, *, plink_prefix=None, genotype_file=None, phenotype_file=None, category_file=None, ld_score_file=None, method='mqs_hew', output_prefix='result', gemma_path=None)`.
- **Input:** A beta/z file, the reference prefix, its PVE estimate and the study genotypes and phenotypes.
- **Output:** `log_values` with the PVE, its standard error and the per-SNP variance.
- **Example code:**

```python
import subprocess
from dnakit.gwas.gemma_p15_cli_functions import *

g, p = 'gwas/geno.txt', 'gwas/pheno.txt'
import tempfile
from pathlib import Path

out = Path(tempfile.mkdtemp())
fit = gemma_summary_variance_components(str(out / "vc"), "gwas/beta.txt", genotype_file=g, phenotype_file=p)
pve = fit["log_values"]["pve_estimates"]
result = gemma_summary_variance_component_errors(str(out / "ci"), "gwas/beta.txt", str(out / "vc" / "result"), pve, genotype_file=g, phenotype_file=p)
print(round(result["log_values"]["pve_estimates"] - pve, 6), result["log_values"]["se_pve"] > 0)
```

- **Example output:**

```text
0.0 True
```

## 6) `dnakit.gwas.gemma_p15_cli_functions.gemma_gxe_association()`

- **Purpose:** Test SNP x environment interaction with a linear mixed model.
- **Method:** Runs `gemma -gxe -lmm`. PLINK files are required because GEMMA 0.98.5 reads BIMBAM genotypes incorrectly in this mode.
- **API:** `dnakit.gwas.gemma_p15_cli_functions.gemma_gxe_association(output_directory, plink_prefix, environment_file, *, kinship_file=None, eigen_values_file=None, eigen_vectors_file=None, phenotype_file=None, covariates_file=None, test='wald', output_prefix='result', miss=None, maf=None, gemma_path=None)`.
- **Input:** PLINK prefix, an environment file, a kinship file (or eigen files), optional phenotype file.
- **Output:** `associations`, `association_count`.
- **Example code:**

```python
import subprocess
from dnakit.gwas.gemma_p15_cli_functions import *

g, p = 'gwas/geno.txt', 'gwas/pheno.txt'
import tempfile
from pathlib import Path

out = Path(tempfile.mkdtemp())
subprocess.run(["gemma", "-g", g, "-p", p, "-gk", "1", "-outdir", str(out), "-o", "kin"], check=True, capture_output=True)
result = gemma_gxe_association(str(out / "gxe"), "gwas/toy", "gwas/env.txt", kinship_file=str(out / "kin.cXX.txt"), phenotype_file=p)
print(result["association_count"], result["associations"][0]["rs"], sorted(result["associations"][0])[:3])
```

- **Example output:**

```text
60 rs1 ['af', 'allele0', 'allele1']
```

## 7) `dnakit.gwas.gemma_p15_cli_functions.gemma_multivariate_association()`

- **Purpose:** Association of SNPs with several phenotypes at once.
- **Method:** Runs `gemma -lmm 1 -n c1 c2 ...` (multivariate linear mixed model, Wald test).
- **API:** `dnakit.gwas.gemma_p15_cli_functions.gemma_multivariate_association(output_directory, phenotype_columns, *, plink_prefix=None, genotype_file=None, phenotype_file=None, annotation_file=None, kinship_file=None, eigen_values_file=None, eigen_vectors_file=None, covariates_file=None, output_prefix='result', miss=None, maf=None, hwe=None, r2=None, gemma_path=None)`.
- **Input:** `phenotype_columns` (at least two), genotypes, phenotypes and a kinship file.
- **Output:** `associations` with one beta per phenotype, `log_values`.
- **Example code:**

```python
import subprocess
from dnakit.gwas.gemma_p15_cli_functions import *

g, p = 'gwas/geno.txt', 'gwas/pheno.txt'
import tempfile
from pathlib import Path

out = Path(tempfile.mkdtemp())
subprocess.run(["gemma", "-g", g, "-p", p, "-gk", "1", "-outdir", str(out), "-o", "kin"], check=True, capture_output=True)
result = gemma_multivariate_association(str(out / "mv"), [1, 2], genotype_file=g, phenotype_file=p, kinship_file=str(out / "kin.cXX.txt"))
print(result["association_count"], [c for c in result["associations"][0] if c.startswith("beta")])
```

- **Example output:**

```text
60 ['beta_1', 'beta_2']
```

## 8) `dnakit.gwas.gemma_p15_cli_functions.gemma_bslmm()`

- **Purpose:** Fit a Bayesian sparse linear mixed model.
- **Method:** Runs `gemma -bslmm 1/2/3` (linear MCMC, ridge/GBLUP, probit). `burn_in` and `sampling_steps` are `-w` and `-s`.
- **API:** `dnakit.gwas.gemma_p15_cli_functions.gemma_bslmm(output_directory, *, model='linear', plink_prefix=None, genotype_file=None, phenotype_file=None, annotation_file=None, phenotype_column=1, kinship_file=None, burn_in=100000, sampling_steps=1000000, max_snps_in_model=None, seed=None, output_prefix='result', miss=None, maf=None, hwe=None, r2=None, gemma_path=None)`.
- **Input:** Genotypes, phenotypes, `model`, MCMC lengths, `seed`.
- **Output:** `snp_effects`, `hyper_parameter_means` (h, pve, rho, pge, pi, n_gamma), files.
- **Example code:**

```python
import subprocess
from dnakit.gwas.gemma_p15_cli_functions import *

g, p = 'gwas/geno.txt', 'gwas/pheno.txt'
import tempfile
from pathlib import Path

out = Path(tempfile.mkdtemp())
result = gemma_bslmm(str(out / "bs"), genotype_file=g, phenotype_file=p, burn_in=200, sampling_steps=1000, seed=1)
print(result["snp_count"], sorted(result["hyper_parameter_means"]), 0 <= result["hyper_parameter_means"]["pve"] <= 1)
ridge = gemma_bslmm(str(out / "ridge"), model="ridge", genotype_file=g, phenotype_file=p)
print(ridge["snp_count"])
```

- **Example output:**

```text
60 ['h', 'n_gamma', 'pge', 'pi', 'pve', 'rho'] True
60
```

## 9) `dnakit.gwas.gemma_p15_cli_functions.gemma_predict_phenotypes()`

- **Purpose:** Predict missing phenotypes from BSLMM estimates or from the kinship matrix.
- **Method:** Runs `gemma -predict` with `-epm -emu [-ebv -k]` (`source='bslmm'`) or with `-k` and `-n` (`source='kinship'`).
- **API:** `dnakit.gwas.gemma_p15_cli_functions.gemma_predict_phenotypes(output_directory, *, source='bslmm', plink_prefix=None, genotype_file=None, phenotype_file=None, parameter_file=None, log_file=None, breeding_value_file=None, kinship_file=None, phenotype_columns=None, probability_scale=False, output_prefix='result', gemma_path=None)`.
- **Input:** The genotype and phenotype files of the fit and the BSLMM parameter and log files, or a kinship file.
- **Output:** `predictions` (one value or `None` per individual), `predicted_count`.
- **Example code:**

```python
import subprocess
from dnakit.gwas.gemma_p15_cli_functions import *

g, p = 'gwas/geno.txt', 'gwas/pheno.txt'
import tempfile
from pathlib import Path

out = Path(tempfile.mkdtemp())
fit = gemma_bslmm(str(out / "bs"), genotype_file=g, phenotype_file=p, burn_in=200, sampling_steps=1000, seed=1)
result = gemma_predict_phenotypes(str(out / "pr"), genotype_file=g, phenotype_file=p,
                                  parameter_file=str(out / "bs" / "result.param.txt"), log_file=str(out / "bs" / "result.log.txt"))
print(result["predicted_count"], result["predictions"][:2], len(result["predictions"]))
```

- **Example output:**

```text
10 [None, None] 80
```

## 10) `dnakit.assembly.quast_p15_cli_functions.quast_assess_assembly()`

- **Purpose:** Evaluate assemblies against a reference: genome fraction, NGAx, misassemblies, mismatches, indels, genes, operons, k-mer metrics, SV comparison and the upper-bound assembly.
- **Method:** Runs `quast.py -r` with minimap2 alignment. `features` is `-g`, `operons` is `--operons`, `k_mer_stats` is `-k`, `sv_bedpe` is `--sv-bedpe`, `upper_bound_assembly` is `--upper-bound-assembly`, `circos` is `--circos`. Plots and HTML are off unless `build_reports`.
- **API:** `dnakit.assembly.quast_p15_cli_functions.quast_assess_assembly(contigs, output_directory, *, reference=None, features=None, operons=None, k_mer_stats=False, k_mer_size=101, upper_bound_assembly=False, upper_bound_reads=None, sv_bedpe=None, circos=False, min_alignment=65, min_identity=95.0, extensive_misassembly_size=1000, eukaryote=False, large=False, min_contig=500, labels=None, build_reports=False, threads=1, quast_path=None)`.
- **Input:** Contig files, a new output folder, `reference` and the options above.
- **Output:** `report`: `{metric: {assembly: value}}` from `report.tsv`.
- **Example code:**

```python
import tempfile
from pathlib import Path

out = Path(tempfile.mkdtemp())
from dnakit.assembly.quast_p15_cli_functions import quast_assess_assembly

result = quast_assess_assembly(["asm/contigs.fa"], str(out / "q"), reference="asm/ref.fa", k_mer_stats=True, k_mer_size=31)
report = result["report"]
for metric in ("Genome fraction (%)", "# misassemblies", "N50", "NGA50", "K-mer-based compl. (%)"):
    print(metric, report[metric]["contigs"])
```

- **Example output:**

```text
Genome fraction (%) 97.5
# misassemblies 0
N50 10500
NGA50 10500
K-mer-based compl. (%) 97.34
```

## 11) `dnakit.assembly.quast_p15_cli_functions.quast_predict_genes()`

- **Purpose:** Predict genes and rRNA genes on contigs.
- **Method:** Runs `quast.py -f` (GeneMarkS, `--mgm` or `--glimmer`) and `--rna-finding` (Barrnap). Needs the gene-finder binaries that QUAST downloads or bundles.
- **API:** `dnakit.assembly.quast_p15_cli_functions.quast_predict_genes(contigs, output_directory, *, gene_finder='genemarks', rna_finding=False, gene_thresholds=(0, 300, 1500, 3000), reference=None, eukaryote=False, min_contig=500, threads=1, quast_path=None)`.
- **Input:** Contig files, a new output folder, `gene_finder`, `rna_finding`.
- **Output:** `report` rows `# predicted genes ...` and `# predicted rRNA genes`.
- **Example code:**

```python
import tempfile
from pathlib import Path

out = Path(tempfile.mkdtemp())
from dnakit.assembly.quast_p15_cli_functions import quast_predict_genes

result = quast_predict_genes(["asm/contigs.fa"], str(out / "q"), gene_finder="mgm", rna_finding=True)
print(result["report"]["# predicted rRNA genes"]["contigs"])
```

- **Example output:**

```text
0 + 0 part
```

## 12) `dnakit.assembly.quast_p15_cli_functions.quast_map_reads()`

- **Purpose:** Map reads to the assembly and reference: coverage, read-based SNPs, indels and structural variants.
- **Method:** Runs QUAST with `--pe1/--pe2`, `--single`, `--pacbio`, `--nanopore`, `--ref-bam` or `--bam`.
- **API:** `dnakit.assembly.quast_p15_cli_functions.quast_map_reads(contigs, output_directory, *, reference=None, pe1=None, pe2=None, single=None, pacbio=None, nanopore=None, reference_bam=None, contig_bam_files=None, detect_variants=True, min_contig=500, threads=1, quast_path=None)`.
- **Input:** Contig files, a new output folder, a reference and reads or alignments.
- **Output:** `report` rows `Mapped (%)`, `Avg. coverage depth`, `# SNPs`.
- **Example code:**

```python
import tempfile
from pathlib import Path

out = Path(tempfile.mkdtemp())
import random
from dnakit.assembly.quast_p15_cli_functions import quast_map_reads

genome = Path("asm/ref.fa").read_text().split("\n")[1]
random.seed(1)
reads = out / "reads.fq"
with reads.open("w") as handle:
    for i in range(2000):
        start = random.randrange(0, len(genome) - 100)
        handle.write(f"@r{i}\n{genome[start:start + 100]}\n+\n{'I' * 100}\n")
result = quast_map_reads(["asm/contigs.fa"], str(out / "q"), reference="asm/ref.fa", single=str(reads))
print(result["report"]["Mapped (%)"]["contigs"] > 90, result["report"]["Avg. coverage depth"]["contigs"])
```

- **Example output:**

```text
True 10
```

## 13) `dnakit.assembly.quast_p15_cli_functions.metaquast_assess_metagenome()`

- **Purpose:** Evaluate metagenome assemblies against several reference genomes, optionally found by BLAST.
- **Method:** Runs `metaquast.py` with `-r` (comma-joined references), `--references-list`, `--max-ref-number` and `--unique-mapping`. The BLAST-based reference search downloads genomes from NCBI and needs network access.
- **API:** `dnakit.assembly.quast_p15_cli_functions.metaquast_assess_metagenome(contigs, output_directory, *, references=None, references_list=None, max_references=50, use_blast_search=True, unique_mapping=False, min_contig=500, threads=1, quast_path=None)`.
- **Input:** Contig files, a new output folder, `references` or `references_list`.
- **Output:** `report` and the `runs_per_reference` folder.
- **Example call (needs the external program; not executed in the documentation build):**

```python
from dnakit.assembly.quast_p15_cli_functions import metaquast_assess_metagenome

result = metaquast_assess_metagenome(["metagenome_contigs.fa"], "metaquast_out", references=["genome_a.fa", "genome_b.fa"], use_blast_search=False)
print(result["report"]["# contigs"])
```

## 14) `dnakit.assembly.busco_p15_cli_functions.busco_assess_completeness()`

- **Purpose:** Assess genome or transcriptome completeness with BUSCO and report assembly statistics and scaffold composition.
- **Method:** Runs `busco -m genome|transcriptome|proteins`. `gene_predictor` selects `--augustus`, `--metaeuk` or `--miniprot` (Prodigal is BUSCO's own choice for prokaryotic lineages); `scaffold_composition` is `--scaffold_composition` and `skip_assembly_statistics` is `--skip_bbtools`.
- **API:** `dnakit.assembly.busco_p15_cli_functions.busco_assess_completeness(input_path, lineage, output_directory, *, mode='genome', gene_predictor=None, output_name='busco', download_path=None, offline=False, auto_lineage=None, evalue=0.001, limit=3, long_mode=False, contig_break=10, scaffold_composition=False, skip_assembly_statistics=False, augustus_species=None, cpus=1, busco_path=None)`.
- **Input:** A FASTA file, a lineage dataset name or folder, a new output folder.
- **Output:** `scores_percent`, `marker_count`, `counts`, `assembly_statistics`, `full_table`, `scaffold_composition`.
- **Example call (needs the external program; not executed in the documentation build):**

```python
from dnakit.assembly.busco_p15_cli_functions import busco_assess_completeness

result = busco_assess_completeness("genome.fa", "bacteria_odb10", "busco_out", gene_predictor="prodigal", scaffold_composition=True, download_path="busco_downloads")
print(result["scores_percent"], result["assembly_statistics"])
```

## 15) `dnakit.variants.delly_p15_cli_functions.delly_genotype_sv_sites()`

- **Purpose:** Genotype known SV sites from short or long reads.
- **Method:** Runs `delly sr -v` or `delly lr -v`.
- **API:** `dnakit.variants.delly_p15_cli_functions.delly_genotype_sv_sites(bam_files, reference, sites_file, output_bcf, *, mode='sr', technology='ont', exclude_regions=None, delly_path=None)`.
- **Input:** Sorted, indexed BAM files, the reference, a site BCF/VCF and a new output BCF.
- **Output:** `variant_count`, `variants` (chromosome, start, end, type, filter, quality).
- **Example code:**

```python
import tempfile
from pathlib import Path

out = Path(tempfile.mkdtemp())
import subprocess
from dnakit.variants.delly_p15_cli_functions import delly_genotype_sv_sites

subprocess.run(["delly", "sr", "-g", "sv/ref.fa", "-o", str(out / "sr.bcf"), "sv/short.bam"], check=True, capture_output=True)
result = delly_genotype_sv_sites(["sv/short.bam"], "sv/ref.fa", str(out / "sr.bcf"), str(out / "genotyped.bcf"))
print(result["variant_count"], result["variants"][0]["type"], result["variants"][0]["start"], result["variants"][0]["end"])
```

- **Example output:**

```text
1 DEL 30000 30801
```

## 16) `dnakit.variants.delly_p15_cli_functions.delly_assembly_sv()`

- **Purpose:** Discover structural variants from assemblies aligned to the reference.
- **Method:** Runs `delly asm` in `squashed` or `diploid` mode.
- **API:** `dnakit.variants.delly_p15_cli_functions.delly_assembly_sv(assembly_bam_files, reference, output_bcf, *, mode='squashed', svtype='ALL', min_mapping_quality=10, min_clip=25, delly_path=None)`.
- **Input:** Assembly BAM files, the reference and a new output BCF.
- **Output:** `variants`, `variant_count`.
- **Example code:**

```python
import tempfile
from pathlib import Path

out = Path(tempfile.mkdtemp())
from dnakit.variants.delly_p15_cli_functions import delly_assembly_sv

result = delly_assembly_sv(["sv/asm.bam"], "sv/ref.fa", str(out / "asm.bcf"))
print(result["variants"][0]["type"], result["variants"][0]["start"], result["variants"][0]["end"])
```

- **Example output:**

```text
DEL 30000 30800
```

## 17) `dnakit.variants.delly_p15_cli_functions.delly_copy_number_variants()`

- **Purpose:** Call copy-number variants from read depth, or re-genotype known CNV sites.
- **Method:** Runs `delly cnv`; with `genotype_sites` it adds `-v`. The coverage and segmentation outputs can feed the read-depth plot.
- **API:** `dnakit.variants.delly_p15_cli_functions.delly_copy_number_variants(bam_file, reference, output_bcf, *, genotype_sites=None, coverage_file=None, segmentation_file=None, exclude_regions=None, min_cnv_size=1000, min_mapping_quality=10, ploidy=2, window=0, delly_path=None)`.
- **Input:** A BAM, the reference, a new BCF and optional coverage and segmentation paths.
- **Output:** `variants`, `regenotyped`, output files.
- **Example code:**

```python
import tempfile
from pathlib import Path

out = Path(tempfile.mkdtemp())
from dnakit.variants.delly_p15_cli_functions import delly_copy_number_variants

result = delly_copy_number_variants("sv/short.bam", "sv/ref.fa", str(out / "cnv.bcf"), coverage_file=str(out / "cov.gz"), segmentation_file=str(out / "seg.bed"))
print(result["status"], result["regenotyped"], Path(out / "cov.gz").is_file())
```

- **Example output:**

```text
completed False True
```

## 18) `dnakit.variants.delly_p15_cli_functions.delly_merge_call_sets()`

- **Purpose:** Merge SV or CNV call sets of several samples.
- **Method:** Runs `delly merge`; `cnv_mode` adds `-e`, `long_reads` adds `-l`.
- **API:** `dnakit.variants.delly_p15_cli_functions.delly_merge_call_sets(bcf_files, output_bcf, *, cnv_mode=False, long_reads=False, reference=None, min_quality=200, min_size=0, max_size=1000000, pass_only=False, delly_path=None)`.
- **Input:** BCF files and a new output BCF.
- **Output:** `variants`, `variant_count`.
- **Example code:**

```python
import tempfile
from pathlib import Path

out = Path(tempfile.mkdtemp())
import subprocess
from dnakit.variants.delly_p15_cli_functions import delly_merge_call_sets

subprocess.run(["delly", "sr", "-g", "sv/ref.fa", "-o", str(out / "a.bcf"), "sv/short.bam"], check=True, capture_output=True)
result = delly_merge_call_sets([str(out / "a.bcf"), str(out / "a.bcf")], str(out / "merged.bcf"), min_quality=0)
print(result["variant_count"])
```

- **Example output:**

```text
1
```

## 19) `dnakit.variants.delly_p15_cli_functions.delly_filter_call_set()`

- **Purpose:** Filter a genotyped SV/CNV call set for somatic or germline variants.
- **Method:** Runs `delly filter -f somatic|germline`; somatic mode needs a sample file with `tumor` and `control` labels.
- **API:** `dnakit.variants.delly_p15_cli_functions.delly_filter_call_set(bcf_file, output_bcf, *, filter_mode='somatic', samples_file=None, min_size=0, max_size=500000000, min_quality=300, min_alt_fraction=0.03, min_genotyped_fraction=0.75, min_coverage=10, max_control_support=0.0, long_reads=False, pass_only=False, tag_only=False, delly_path=None)`.
- **Input:** A BCF, a new output BCF, `filter_mode` and thresholds.
- **Output:** `variants`, `variant_count`.
- **Example code:**

```python
import tempfile
from pathlib import Path

out = Path(tempfile.mkdtemp())
import subprocess
from dnakit.variants.delly_p15_cli_functions import delly_filter_call_set

subprocess.run(["delly", "sr", "-g", "sv/ref.fa", "-o", str(out / "a.bcf"), "sv/short.bam"], check=True, capture_output=True)
result = delly_filter_call_set(str(out / "a.bcf"), str(out / "filtered.bcf"), filter_mode="germline", min_quality=0)
print(result["status"], result["filter_mode"])
```

- **Example output:**

```text
completed germline
```

## 20) `dnakit.variants.delly_p15_cli_functions.delly_translocations_to_breakends()`

- **Purpose:** Convert DELLY translocation records into pairs of breakend (BND) records.
- **Method:** Runs DELLY's `scripts/delly2bnd.py` (the file from the DELLY repository, passed as `script_path`) with a Python that has cyvcf2.
- **API:** `dnakit.variants.delly_p15_cli_functions.delly_translocations_to_breakends(vcf_file, reference, output_vcf, script_path, *, python_path=None)`.
- **Input:** A DELLY VCF/BCF, the reference, a new VCF and the script path.
- **Output:** `breakend_record_count`.
- **Example call (needs the external program; not executed in the documentation build):**

```python
from dnakit.variants.delly_p15_cli_functions import delly_translocations_to_breakends

result = delly_translocations_to_breakends("delly_calls.bcf", "ref.fa", "breakends.vcf", "delly/scripts/delly2bnd.py")
print(result["breakend_record_count"])
```

## 21) `dnakit.variants.delly_p15_cli_functions.delly_plot_read_depth()`

- **Purpose:** Plot the read-depth copy-number profile and its segmentation.
- **Method:** Runs DELLY's `R/rd.R` with `Rscript` (packages DNAcopy, ggplot2, scales, gtable); without a segmentation file it segments by circular binary segmentation.
- **API:** `dnakit.variants.delly_p15_cli_functions.delly_plot_read_depth(coverage_file, output_directory, script_path, *, segmentation_file=None, rscript_path=None)`.
- **Input:** The coverage file of `delly_copy_number_variants`, the script path and a new output folder.
- **Output:** `images` (PNG paths) and `segmentation`.
- **Example call (needs the external program; not executed in the documentation build):**

```python
from dnakit.variants.delly_p15_cli_functions import delly_plot_read_depth

result = delly_plot_read_depth("cov.gz", "rd_plots", "delly/R/rd.R")
print(result["images"])
```
