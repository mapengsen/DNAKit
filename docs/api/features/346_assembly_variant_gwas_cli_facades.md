# 组装、变异与 GWAS 命令行适配器：GEMMA、QUAST、BUSCO 与 DELLY

这一组函数通过固定的参数白名单、不经过 shell 调用四个由使用者自行安装的程序：[GEMMA](https://github.com/genetics-statistics/GEMMA)（由汇总统计估计方差组分、特征分解、G x E、多变量 LMM、BSLMM、预测）、[QUAST](https://github.com/ablab/quast)（基于参考的组装评估、基因与 rRNA 预测、reads 比对、MetaQUAST）、[BUSCO](https://busco.ezlab.org/)（指定基因预测器的完整性评估、组装统计、scaffold 组成）和 [DELLY](https://github.com/dellytools/delly)（基因分型、基于组装的 SV、拷贝数、合并、过滤、断点格式转换、读深图）。缺少程序时报错并给出安装命令；输出文件和文件夹必须尚不存在。

测试版本：GEMMA 0.98.5、QUAST 5.2.0、DELLY 2.7.0。注意：GEMMA 0.98.5 在 `-gs` 写完文件后会因堆错误中止（适配器接受这种情况并在 `warnings` 中说明）；其 `-gw`（LD 分数权重）上游尚未实现，BIMBAM 的 `-gxe` 会读错文件，所以 G x E 函数使用 PLINK 文件。测试环境没有 BUSCO 和 R：BUSCO 适配器用模拟其命令行和输出文件的替身程序测试，读深图（`R/rd.R`）有文档但未运行。因此 BUSCO 与读深图的示例调用不执行。

- **示例代码：**

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

- **示例输出：**

```text
['asm.bam', 'long.bam', 'short.bam'] 7 ['contigs.fa', 'ref.fa']
```

## 1) `dnakit.gwas.gemma_p15_cli_functions.gemma_snp_covariance_matrix()`

- **作用：** 为基于汇总统计的方差组分计算 SNP 方差-协方差（S）矩阵。
- **处理方法：** 运行 `gemma -gs`。用有表型的个体估计每个类别 SNP 的 S 矩阵。GEMMA 0.98.5 写完文件后会中止；这被接受并在 `warnings` 中说明。
- **API：** `dnakit.gwas.gemma_p15_cli_functions.gemma_snp_covariance_matrix(output_directory, *, plink_prefix=None, genotype_file=None, phenotype_file=None, category_file=None, sample=None, output_prefix='result', miss=None, maf=None, hwe=None, r2=None, not_snp=False, gemma_path=None)`。
- **输入：** BIMBAM 或 PLINK 基因型、表型文件，可选 `category_file`、`sample`。
- **输出：** `s_matrix`、`sizes`、`warnings`，以及新文件夹中的文件。
- **示例代码：**

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

- **示例输出：**

```text
[60.0, 70.0] completed 1
0.01604
```

## 2) `dnakit.gwas.gemma_p15_cli_functions.gemma_summary_statistic_q_vector()`

- **作用：** 由 GWAS 汇总统计文件计算 q 向量。
- **处理方法：** 运行 `gemma -gq`。GEMMA 0.98.5 还需要一个表型文件才能启动，其数值不参与计算。
- **API：** `dnakit.gwas.gemma_p15_cli_functions.gemma_summary_statistic_q_vector(output_directory, beta_file, phenotype_file, *, category_file=None, output_prefix='result', gemma_path=None)`。
- **输入：** beta/z 文件（SNP、N、Z、等位基因）、表型文件、可选类别文件。
- **输出：** `q`、`q_errors`、`sizes`。
- **示例代码：**

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

- **示例输出：**

```text
[55.0, 70.0] [-0.00796]
```

## 3) `dnakit.gwas.gemma_p15_cli_functions.gemma_eigen_decomposition()`

- **作用：** 对亲缘（相关性）矩阵做特征分解。
- **处理方法：** 运行 `gemma -eigen`：取被分析个体对应的行列，中心化后分解。
- **API：** `dnakit.gwas.gemma_p15_cli_functions.gemma_eigen_decomposition(output_directory, kinship_file, *, plink_prefix=None, genotype_file=None, phenotype_file=None, covariates_file=None, output_prefix='result', gemma_path=None)`。
- **输入：** 由 `gemma -gk` 得到的亲缘矩阵、基因型和表型文件、新的输出文件夹。
- **输出：** `eigen_values`、`individual_count`；`<prefix>.eigenU.txt` 为特征向量。
- **示例代码：**

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

- **示例输出：**

```text
70 1.7967
```

## 4) `dnakit.gwas.gemma_p15_cli_functions.gemma_summary_variance_components()`

- **作用：** 由 GWAS 汇总统计和参考 LD 面板估计方差组分与遗传度。
- **处理方法：** 运行 `gemma -vc`，算法为 MQS-HEW（`-vc 1`）或 MQS-LDW（`-vc 2`，需要 LD 分数文件）。
- **API：** `dnakit.gwas.gemma_p15_cli_functions.gemma_summary_variance_components(output_directory, beta_file, *, plink_prefix=None, genotype_file=None, phenotype_file=None, category_file=None, ld_score_file=None, method='mqs_hew', covariates_file=None, sample=None, output_prefix='result', miss=None, maf=None, hwe=None, r2=None, not_snp=False, gemma_path=None)`。
- **输入：** beta/z 文件、基因型（PLINK 或 BIMBAM）、表型文件，可选类别与 `method`。
- **输出：** `log_values`，含 `pve_estimates`、`se_pve`、`sigma2_estimates`、`enrichment`。
- **示例代码：**

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

- **示例输出：**

```text
55.0 -0.4961 0.0203
```

## 5) `dnakit.gwas.gemma_p15_cli_functions.gemma_summary_variance_component_errors()`

- **作用：** 用渐近方法和参考面板计算汇总统计方差组分的标准误。
- **处理方法：** 运行带 `-ref`（先前拟合的前缀）和 `-pve` 的 `gemma -ci`。
- **API：** `dnakit.gwas.gemma_p15_cli_functions.gemma_summary_variance_component_errors(output_directory, beta_file, reference_prefix, pve, *, plink_prefix=None, genotype_file=None, phenotype_file=None, category_file=None, ld_score_file=None, method='mqs_hew', output_prefix='result', gemma_path=None)`。
- **输入：** beta/z 文件、参考前缀、其 PVE 估计值以及研究的基因型和表型。
- **输出：** `log_values`，含 PVE、其标准误以及每个 SNP 的方差。
- **示例代码：**

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

- **示例输出：**

```text
0.0 True
```

## 6) `dnakit.gwas.gemma_p15_cli_functions.gemma_gxe_association()`

- **作用：** 用线性混合模型检验 SNP 与环境的交互作用。
- **处理方法：** 运行 `gemma -gxe -lmm`。必须使用 PLINK 文件，因为 GEMMA 0.98.5 在此模式下会错误读取 BIMBAM 基因型。
- **API：** `dnakit.gwas.gemma_p15_cli_functions.gemma_gxe_association(output_directory, plink_prefix, environment_file, *, kinship_file=None, eigen_values_file=None, eigen_vectors_file=None, phenotype_file=None, covariates_file=None, test='wald', output_prefix='result', miss=None, maf=None, gemma_path=None)`。
- **输入：** PLINK 前缀、环境文件、亲缘矩阵（或特征分解文件），可选表型文件。
- **输出：** `associations`、`association_count`。
- **示例代码：**

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

- **示例输出：**

```text
60 rs1 ['af', 'allele0', 'allele1']
```

## 7) `dnakit.gwas.gemma_p15_cli_functions.gemma_multivariate_association()`

- **作用：** 同时检验 SNP 与多个表型的关联。
- **处理方法：** 运行 `gemma -lmm 1 -n c1 c2 ...`（多变量线性混合模型，Wald 检验）。
- **API：** `dnakit.gwas.gemma_p15_cli_functions.gemma_multivariate_association(output_directory, phenotype_columns, *, plink_prefix=None, genotype_file=None, phenotype_file=None, annotation_file=None, kinship_file=None, eigen_values_file=None, eigen_vectors_file=None, covariates_file=None, output_prefix='result', miss=None, maf=None, hwe=None, r2=None, gemma_path=None)`。
- **输入：** `phenotype_columns`（至少两个）、基因型、表型和亲缘矩阵。
- **输出：** `associations`（每个表型一个 beta）、`log_values`。
- **示例代码：**

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

- **示例输出：**

```text
60 ['beta_1', 'beta_2']
```

## 8) `dnakit.gwas.gemma_p15_cli_functions.gemma_bslmm()`

- **作用：** 拟合贝叶斯稀疏线性混合模型。
- **处理方法：** 运行 `gemma -bslmm 1/2/3`（线性 MCMC、岭回归/GBLUP、probit）。`burn_in` 与 `sampling_steps` 即 `-w` 和 `-s`。
- **API：** `dnakit.gwas.gemma_p15_cli_functions.gemma_bslmm(output_directory, *, model='linear', plink_prefix=None, genotype_file=None, phenotype_file=None, annotation_file=None, phenotype_column=1, kinship_file=None, burn_in=100000, sampling_steps=1000000, max_snps_in_model=None, seed=None, output_prefix='result', miss=None, maf=None, hwe=None, r2=None, gemma_path=None)`。
- **输入：** 基因型、表型、`model`、MCMC 长度、`seed`。
- **输出：** `snp_effects`、`hyper_parameter_means`（h、pve、rho、pge、pi、n_gamma）、文件。
- **示例代码：**

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

- **示例输出：**

```text
60 ['h', 'n_gamma', 'pge', 'pi', 'pve', 'rho'] True
60
```

## 9) `dnakit.gwas.gemma_p15_cli_functions.gemma_predict_phenotypes()`

- **作用：** 由 BSLMM 估计值或亲缘矩阵预测缺失表型。
- **处理方法：** 运行 `gemma -predict`，参数为 `-epm -emu [-ebv -k]`（`source='bslmm'`）或 `-k` 与 `-n`（`source='kinship'`）。
- **API：** `dnakit.gwas.gemma_p15_cli_functions.gemma_predict_phenotypes(output_directory, *, source='bslmm', plink_prefix=None, genotype_file=None, phenotype_file=None, parameter_file=None, log_file=None, breeding_value_file=None, kinship_file=None, phenotype_columns=None, probability_scale=False, output_prefix='result', gemma_path=None)`。
- **输入：** 拟合所用的基因型与表型文件及 BSLMM 的参数和日志文件，或亲缘矩阵。
- **输出：** `predictions`（每个个体一个值或 `None`）、`predicted_count`。
- **示例代码：**

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

- **示例输出：**

```text
10 [None, None] 80
```

## 10) `dnakit.assembly.quast_p15_cli_functions.quast_assess_assembly()`

- **作用：** 对照参考序列评估组装：基因组覆盖率、NGAx、错误组装、错配、插入缺失、基因、操纵子、k-mer 指标、SV 比较和上限组装。
- **处理方法：** 运行 `quast.py -r`，用 minimap2 比对。`features` 即 `-g`，`operons` 即 `--operons`，`k_mer_stats` 即 `-k`，`sv_bedpe` 即 `--sv-bedpe`，`upper_bound_assembly` 即 `--upper-bound-assembly`，`circos` 即 `--circos`。除非 `build_reports`，不生成图和 HTML。
- **API：** `dnakit.assembly.quast_p15_cli_functions.quast_assess_assembly(contigs, output_directory, *, reference=None, features=None, operons=None, k_mer_stats=False, k_mer_size=101, upper_bound_assembly=False, upper_bound_reads=None, sv_bedpe=None, circos=False, min_alignment=65, min_identity=95.0, extensive_misassembly_size=1000, eukaryote=False, large=False, min_contig=500, labels=None, build_reports=False, threads=1, quast_path=None)`。
- **输入：** contig 文件、新的输出文件夹、`reference` 以及上述选项。
- **输出：** `report`：来自 `report.tsv` 的 `{指标: {组装: 值}}`。
- **示例代码：**

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

- **示例输出：**

```text
Genome fraction (%) 97.5
# misassemblies 0
N50 10500
NGA50 10500
K-mer-based compl. (%) 97.34
```

## 11) `dnakit.assembly.quast_p15_cli_functions.quast_predict_genes()`

- **作用：** 在 contig 上预测基因与 rRNA 基因。
- **处理方法：** 运行 `quast.py -f`（GeneMarkS、`--mgm` 或 `--glimmer`）和 `--rna-finding`（Barrnap）。需要 QUAST 自带或下载的基因预测程序。
- **API：** `dnakit.assembly.quast_p15_cli_functions.quast_predict_genes(contigs, output_directory, *, gene_finder='genemarks', rna_finding=False, gene_thresholds=(0, 300, 1500, 3000), reference=None, eukaryote=False, min_contig=500, threads=1, quast_path=None)`。
- **输入：** contig 文件、新的输出文件夹、`gene_finder`、`rna_finding`。
- **输出：** `report` 中的 `# predicted genes ...` 与 `# predicted rRNA genes` 行。
- **示例代码：**

```python
import tempfile
from pathlib import Path

out = Path(tempfile.mkdtemp())
from dnakit.assembly.quast_p15_cli_functions import quast_predict_genes

result = quast_predict_genes(["asm/contigs.fa"], str(out / "q"), gene_finder="mgm", rna_finding=True)
print(result["report"]["# predicted rRNA genes"]["contigs"])
```

- **示例输出：**

```text
0 + 0 part
```

## 12) `dnakit.assembly.quast_p15_cli_functions.quast_map_reads()`

- **作用：** 把 reads 比对到组装和参考序列：覆盖度、基于 reads 的 SNP、插入缺失与结构变异。
- **处理方法：** 运行带 `--pe1/--pe2`、`--single`、`--pacbio`、`--nanopore`、`--ref-bam` 或 `--bam` 的 QUAST。
- **API：** `dnakit.assembly.quast_p15_cli_functions.quast_map_reads(contigs, output_directory, *, reference=None, pe1=None, pe2=None, single=None, pacbio=None, nanopore=None, reference_bam=None, contig_bam_files=None, detect_variants=True, min_contig=500, threads=1, quast_path=None)`。
- **输入：** contig 文件、新的输出文件夹、参考序列以及 reads 或比对文件。
- **输出：** `report` 中的 `Mapped (%)`、`Avg. coverage depth`、`# SNPs` 行。
- **示例代码：**

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

- **示例输出：**

```text
True 10
```

## 13) `dnakit.assembly.quast_p15_cli_functions.metaquast_assess_metagenome()`

- **作用：** 对照多个参考基因组评估宏基因组组装，参考基因组可由 BLAST 搜索得到。
- **处理方法：** 运行 `metaquast.py`，参数为 `-r`（逗号连接的参考序列）、`--references-list`、`--max-ref-number` 与 `--unique-mapping`。基于 BLAST 的参考序列搜索会从 NCBI 下载基因组，需要网络。
- **API：** `dnakit.assembly.quast_p15_cli_functions.metaquast_assess_metagenome(contigs, output_directory, *, references=None, references_list=None, max_references=50, use_blast_search=True, unique_mapping=False, min_contig=500, threads=1, quast_path=None)`。
- **输入：** contig 文件、新的输出文件夹、`references` 或 `references_list`。
- **输出：** `report` 与 `runs_per_reference` 文件夹。
- **示例调用（需要外部程序，文档构建时不执行）：**

```python
from dnakit.assembly.quast_p15_cli_functions import metaquast_assess_metagenome

result = metaquast_assess_metagenome(["metagenome_contigs.fa"], "metaquast_out", references=["genome_a.fa", "genome_b.fa"], use_blast_search=False)
print(result["report"]["# contigs"])
```

## 14) `dnakit.assembly.busco_p15_cli_functions.busco_assess_completeness()`

- **作用：** 用 BUSCO 评估基因组或转录组完整性，并报告组装统计和 scaffold 组成。
- **处理方法：** 运行 `busco -m genome|transcriptome|proteins`。`gene_predictor` 选择 `--augustus`、`--metaeuk` 或 `--miniprot`（原核谱系 BUSCO 自动使用 Prodigal）；`scaffold_composition` 即 `--scaffold_composition`，`skip_assembly_statistics` 即 `--skip_bbtools`。
- **API：** `dnakit.assembly.busco_p15_cli_functions.busco_assess_completeness(input_path, lineage, output_directory, *, mode='genome', gene_predictor=None, output_name='busco', download_path=None, offline=False, auto_lineage=None, evalue=0.001, limit=3, long_mode=False, contig_break=10, scaffold_composition=False, skip_assembly_statistics=False, augustus_species=None, cpus=1, busco_path=None)`。
- **输入：** FASTA 文件、谱系数据集名称或文件夹、新的输出文件夹。
- **输出：** `scores_percent`、`marker_count`、`counts`、`assembly_statistics`、`full_table`、`scaffold_composition`。
- **示例调用（需要外部程序，文档构建时不执行）：**

```python
from dnakit.assembly.busco_p15_cli_functions import busco_assess_completeness

result = busco_assess_completeness("genome.fa", "bacteria_odb10", "busco_out", gene_predictor="prodigal", scaffold_composition=True, download_path="busco_downloads")
print(result["scores_percent"], result["assembly_statistics"])
```

## 15) `dnakit.variants.delly_p15_cli_functions.delly_genotype_sv_sites()`

- **作用：** 用短读长或长读长对已知 SV 位点分型。
- **处理方法：** 运行 `delly sr -v` 或 `delly lr -v`。
- **API：** `dnakit.variants.delly_p15_cli_functions.delly_genotype_sv_sites(bam_files, reference, sites_file, output_bcf, *, mode='sr', technology='ont', exclude_regions=None, delly_path=None)`。
- **输入：** 已排序并建索引的 BAM、参考序列、位点 BCF/VCF 和新的输出 BCF。
- **输出：** `variant_count`、`variants`（染色体、起止、类型、过滤、质量）。
- **示例代码：**

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

- **示例输出：**

```text
1 DEL 30000 30801
```

## 16) `dnakit.variants.delly_p15_cli_functions.delly_assembly_sv()`

- **作用：** 由比对到参考序列的组装发现结构变异。
- **处理方法：** 以 `squashed` 或 `diploid` 模式运行 `delly asm`。
- **API：** `dnakit.variants.delly_p15_cli_functions.delly_assembly_sv(assembly_bam_files, reference, output_bcf, *, mode='squashed', svtype='ALL', min_mapping_quality=10, min_clip=25, delly_path=None)`。
- **输入：** 组装 BAM、参考序列和新的输出 BCF。
- **输出：** `variants`、`variant_count`。
- **示例代码：**

```python
import tempfile
from pathlib import Path

out = Path(tempfile.mkdtemp())
from dnakit.variants.delly_p15_cli_functions import delly_assembly_sv

result = delly_assembly_sv(["sv/asm.bam"], "sv/ref.fa", str(out / "asm.bcf"))
print(result["variants"][0]["type"], result["variants"][0]["start"], result["variants"][0]["end"])
```

- **示例输出：**

```text
DEL 30000 30800
```

## 17) `dnakit.variants.delly_p15_cli_functions.delly_copy_number_variants()`

- **作用：** 由读深调用拷贝数变异，或对已知 CNV 位点重新分型。
- **处理方法：** 运行 `delly cnv`；给出 `genotype_sites` 时加上 `-v`。覆盖度与分段输出可用于读深图。
- **API：** `dnakit.variants.delly_p15_cli_functions.delly_copy_number_variants(bam_file, reference, output_bcf, *, genotype_sites=None, coverage_file=None, segmentation_file=None, exclude_regions=None, min_cnv_size=1000, min_mapping_quality=10, ploidy=2, window=0, delly_path=None)`。
- **输入：** BAM、参考序列、新的 BCF，以及可选的覆盖度和分段输出路径。
- **输出：** `variants`、`regenotyped`、输出文件。
- **示例代码：**

```python
import tempfile
from pathlib import Path

out = Path(tempfile.mkdtemp())
from dnakit.variants.delly_p15_cli_functions import delly_copy_number_variants

result = delly_copy_number_variants("sv/short.bam", "sv/ref.fa", str(out / "cnv.bcf"), coverage_file=str(out / "cov.gz"), segmentation_file=str(out / "seg.bed"))
print(result["status"], result["regenotyped"], Path(out / "cov.gz").is_file())
```

- **示例输出：**

```text
completed False True
```

## 18) `dnakit.variants.delly_p15_cli_functions.delly_merge_call_sets()`

- **作用：** 合并多个样本的 SV 或 CNV 结果。
- **处理方法：** 运行 `delly merge`；`cnv_mode` 加 `-e`，`long_reads` 加 `-l`。
- **API：** `dnakit.variants.delly_p15_cli_functions.delly_merge_call_sets(bcf_files, output_bcf, *, cnv_mode=False, long_reads=False, reference=None, min_quality=200, min_size=0, max_size=1000000, pass_only=False, delly_path=None)`。
- **输入：** BCF 文件和新的输出 BCF。
- **输出：** `variants`、`variant_count`。
- **示例代码：**

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

- **示例输出：**

```text
1
```

## 19) `dnakit.variants.delly_p15_cli_functions.delly_filter_call_set()`

- **作用：** 对已分型的 SV/CNV 结果做体细胞或胚系过滤。
- **处理方法：** 运行 `delly filter -f somatic|germline`；体细胞模式需要标注 `tumor` 与 `control` 的样本文件。
- **API：** `dnakit.variants.delly_p15_cli_functions.delly_filter_call_set(bcf_file, output_bcf, *, filter_mode='somatic', samples_file=None, min_size=0, max_size=500000000, min_quality=300, min_alt_fraction=0.03, min_genotyped_fraction=0.75, min_coverage=10, max_control_support=0.0, long_reads=False, pass_only=False, tag_only=False, delly_path=None)`。
- **输入：** BCF、新的输出 BCF、`filter_mode` 与阈值。
- **输出：** `variants`、`variant_count`。
- **示例代码：**

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

- **示例输出：**

```text
completed germline
```

## 20) `dnakit.variants.delly_p15_cli_functions.delly_translocations_to_breakends()`

- **作用：** 把 DELLY 的易位记录转换为成对的断点（BND）记录。
- **处理方法：** 用带 cyvcf2 的 Python 运行 DELLY 的 `scripts/delly2bnd.py`（来自 DELLY 仓库，通过 `script_path` 传入）。
- **API：** `dnakit.variants.delly_p15_cli_functions.delly_translocations_to_breakends(vcf_file, reference, output_vcf, script_path, *, python_path=None)`。
- **输入：** DELLY 的 VCF/BCF、参考序列、新的 VCF 和脚本路径。
- **输出：** `breakend_record_count`。
- **示例调用（需要外部程序，文档构建时不执行）：**

```python
from dnakit.variants.delly_p15_cli_functions import delly_translocations_to_breakends

result = delly_translocations_to_breakends("delly_calls.bcf", "ref.fa", "breakends.vcf", "delly/scripts/delly2bnd.py")
print(result["breakend_record_count"])
```

## 21) `dnakit.variants.delly_p15_cli_functions.delly_plot_read_depth()`

- **作用：** 绘制读深拷贝数曲线及其分段。
- **处理方法：** 用 `Rscript` 运行 DELLY 的 `R/rd.R`（需要 DNAcopy、ggplot2、scales、gtable 包）；没有分段文件时用循环二元分割进行分段。
- **API：** `dnakit.variants.delly_p15_cli_functions.delly_plot_read_depth(coverage_file, output_directory, script_path, *, segmentation_file=None, rscript_path=None)`。
- **输入：** `delly_copy_number_variants` 输出的覆盖度文件、脚本路径和新的输出文件夹。
- **输出：** `images`（PNG 路径）与 `segmentation`。
- **示例调用（需要外部程序，文档构建时不执行）：**

```python
from dnakit.variants.delly_p15_cli_functions import delly_plot_read_depth

result = delly_plot_read_depth("cov.gz", "rd_plots", "delly/R/rd.R")
print(result["images"])
```
