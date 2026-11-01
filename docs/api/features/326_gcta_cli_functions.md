# GCTA：关系矩阵、REML、关联分析与汇总统计分析

这一组函数是 `dnakit.gwas` 中针对 [GCTA](https://yanglab.westlake.edu.cn/software/gcta/)（Genome-wide Complex Trait Analysis，Yang Lab；以 MIT 许可发布二进制程序）的命令行封装。每个函数用明确的参数列表（不经过 shell）运行一个固定的 GCTA 分析，拒绝覆盖已存在的结果，并返回包含 `status`、`backend`、`version`、`command`、`return_code`、`outputs`（结果名到文件）、`log_tail` 和主结果文件小结 `summary` 的字典。只有退出码为 0 且所有预期输出都存在时 `status` 才是 `completed`。值为 `None` 的参数不会传给 GCTA，所以使用 GCTA 自己的默认值。基因型为 PLINK 二进制文件集。测试版本为 GCTA 1.95.4；fastGWA-GE 的 `envir_file` 需要 1.95 或更高版本，COJO 结果文件的 `.cojo` 后缀在 1.94 与 1.95 之间不同（两者都能识别）。

GCTA 是外部程序。请先安装（发布的二进制文件，或用 `conda install -c bioconda gcta` 安装 1.94.1），并把 `gcta64` 放到 PATH 或传入 `gcta_path`；找不到时函数抛出带安装提示的 `BackendUnavailableError`。`gcta_simulate_phenotypes` 每次运行都重新抽随机数，所以每次调用得到的表型不同。示例直接从模块导入，所有函数都以 `verbose=False` 调用。

## 示例数据

示例读取下面这些很小的文件。函数从不覆盖已有文件，所以每个输出都写入新建的临时目录。

唯一的输入由第一步生成：由固定种子的随机数生成的 PLINK 文件集 `toy`（200 个样本、2600 个 SNP，位于 1、2 号和 X 染色体）、含 60 个因果 SNP 的数量性状、与之相关的第二个性状、二分类性状、协变量和环境文件，以及很小的 SNP、基因和集合列表。

- **示例代码：**

```python
import numpy as np
from pathlib import Path

rng = np.random.default_rng(1)
samples, per_chromosome = 200, (1300, 1200, 100)
snps = sum(per_chromosome)
frequency = rng.uniform(0.1, 0.5, snps)
genotypes = rng.binomial(2, frequency, size=(samples, snps))
chromosomes = [1] * per_chromosome[0] + [2] * per_chromosome[1] + [23] * per_chromosome[2]
positions = [1000 * (k + 1) for count in per_chromosome for k in range(count)]
code = {0: 0b00, 1: 0b10, 2: 0b11}
with open("toy.bed", "wb") as handle:
    handle.write(bytes([0x6C, 0x1B, 0x01]))
    for j in range(snps):
        packed = bytearray((samples + 3) // 4)
        for i, value in enumerate(genotypes[:, j]):
            packed[i // 4] |= code[int(value)] << (2 * (i % 4))
        handle.write(bytes(packed))
Path("toy.bim").write_text("".join(f"{chromosomes[j]}\trs{j + 1}\t0\t{positions[j]}\tA\tG\n" for j in range(snps)))
Path("toy.fam").write_text("".join(f"F{i + 1} I{i + 1} 0 0 {1 + i % 2} -9\n" for i in range(samples)))
beta = np.zeros(snps)
causal = rng.choice(snps, 60, replace=False)
beta[causal] = rng.normal(0, 0.3, 60)
standardised = (genotypes - 2 * frequency) / np.sqrt(2 * frequency * (1 - frequency))
trait = standardised @ beta + rng.normal(0, 1.0, samples)
second = 0.5 * trait + rng.normal(0, 1.0, samples)
binary = (trait > np.median(trait)).astype(int)
covariate = rng.normal(0, 1.0, samples)
lines = {
    "toy.phen": [f"F{i + 1} I{i + 1} {trait[i]:.5f}" for i in range(samples)],
    "toy2.phen": [f"F{i + 1} I{i + 1} {trait[i]:.5f} {second[i]:.5f}" for i in range(samples)],
    "toy.bin.phen": [f"F{i + 1} I{i + 1} {binary[i]}" for i in range(samples)],
    "toy.qcovar": [f"F{i + 1} I{i + 1} {covariate[i]:.4f}" for i in range(samples)],
    "toy.envir": [f"F{i + 1} I{i + 1} {covariate[(i + 7) % samples]:.4f}" for i in range(samples)],
    "toy.env_class": [f"F{i + 1} I{i + 1} {i % 3}" for i in range(samples)],
    "toy.subpop": [f"F{i + 1} I{i + 1} {'PopA' if i % 2 == 0 else 'PopB'}" for i in range(samples)],
    "keep150.txt": [f"F{i + 1} I{i + 1}" for i in range(150)],
    "causal.snplist": ["rs5 0.1", "rs100 -0.2"],
    "ld.snplist": ["rs5", "rs100"],
    "joint.snplist": ["rs5", "rs100", "rs300"],
    "cond.snplist": ["rs5"],
    "genes.txt": ["1 1 400000 Gene1", "1 400001 800000 Gene2", "2 1 400000 Gene3"],
    "genes_burden.txt": ["gene chr start end strand", "G1 1 1 100000 +", "G2 1 200000 400000 +", "G3 2 1 300000 -"],
    "genes_acat.txt": ["1 1000 100000 G1", "1 200000 400000 G2", "2 1000 300000 G3"],
    "sets.txt": ["S1", "rs1", "rs2", "rs3", "rs4", "END", "S2", "rs10", "rs11", "rs12", "END"],
    "set1.snps": [f"rs{k}" for k in range(1, 301)],
    "set2.snps": [f"rs{k}" for k in range(301, 701)],
}
for name, content in lines.items():
    Path(name).write_text("\n".join(content) + "\n")
print(snps, samples, len(lines))
```

- **示例输出：**

```text
2600 200 18
```

## 1) `dnakit.gwas.gcta_recode_standardized_genotypes()`

- **作用：** 把 PLINK 基因型标准化后以压缩文本导出。
- **处理方法：** 运行 `gcta64 --recode-std`，每个基因型标准化为 (x-2p)/sqrt(2p(1-p))。
- **API：** `dnakit.gwas.gcta_recode_standardized_genotypes(bfile, out_prefix, *, keep=None, extract=None, chromosome=None, maf=None, max_maf=None, threads=None, ..., verbose=True)`。
- **输入：** PLINK 前缀和新的输出前缀。
- **输出：** 包含 `status`、`backend`、`version`、`command`、`return_code`、`outputs`、`log_tail` 的字典，另有`n_individuals`、`n_snps`；矩阵在 `<前缀>.xmat.gz`。
- **示例代码：**

```python
import tempfile
from pathlib import Path

out = Path(tempfile.mkdtemp())
from dnakit.gwas.gcta_p08_cli_functions import gcta_recode_standardized_genotypes
import gzip

result = gcta_recode_standardized_genotypes("toy", str(out / "rs"), chromosome=1, verbose=False)
with gzip.open(result["outputs"]["xmat"], "rt") as handle:
    rows = [line.split() for line in handle]
print(result["status"], result["n_individuals"], result["n_snps"], len(rows) - 2)
```

- **示例输出：**

```text
completed 200 1300 200
```

## 2) `dnakit.gwas.gcta_make_grm()`

- **作用：** 由 SNP 基因型计算稠密遗传关系矩阵（GRM）。
- **处理方法：** `grm_type` 选择 `--make-grm`（常染色体加性）、`--make-grm-xchr`、`--make-grm-d`（显性）或 `--make-grm-inbred`。
- **API：** `dnakit.gwas.gcta_make_grm(out_prefix, *, bfile=None, mbfile=None, grm_type='autosome', keep=None, extract=None, chromosome=None, ..., verbose=True)`。
- **输入：** PLINK 前缀（`bfile`）或前缀列表文件（`mbfile`）、新的输出前缀和 `grm_type`。
- **输出：** 包含 `status`、`backend`、`version`、`command`、`return_code`、`outputs`、`log_tail` 的字典，另有`grm_prefix`、`n_individuals`；文件 `.grm.bin`、`.grm.N.bin`、`.grm.id`。
- **示例代码：**

```python
import tempfile
from pathlib import Path

out = Path(tempfile.mkdtemp())
from dnakit.gwas.gcta_p08_cli_functions import gcta_make_grm
import struct

result = gcta_make_grm(str(out / "g"), bfile="toy", verbose=False)
raw = Path(result["grm_prefix"] + ".grm.bin").read_bytes()
values = struct.unpack(f"<{len(raw) // 4}f", raw)
diagonal = [values[i * (i + 1) // 2 + i] for i in range(200)]
print(result["status"], result["n_individuals"], len(values), round(sum(diagonal) / 200, 2))
dominance = gcta_make_grm(str(out / "d"), bfile="toy", grm_type="dominance", verbose=False)
print(Path(dominance["grm_prefix"]).name, sorted(Path(p).name for p in dominance["outputs"].values() if "grm" in p))
```

- **示例输出：**

```text
completed 200 20100 1.0
d.d ['d.d.grm.N.bin', 'd.d.grm.bin', 'd.d.grm.id']
```

## 3) `dnakit.gwas.gcta_make_sparse_grm()`

- **作用：** 把稠密 GRM 按阈值转换成 fastGWA 使用的稀疏 GRM。
- **处理方法：** 运行 `--make-bK-sparse`（稀疏格式）或 `--make-bK`（阈值化后的稠密格式），小于 `cutoff` 的非对角元素置零。
- **API：** `dnakit.gwas.gcta_make_sparse_grm(grm_prefix, out_prefix, *, cutoff=0.05, output_format='sparse', gcta_path=None, timeout_seconds=3600.0, verbose=True)`。
- **输入：** GRM 前缀、新的输出前缀和 `cutoff`（默认 0.05）。
- **输出：** 包含 `status`、`backend`、`version`、`command`、`return_code`、`outputs`、`log_tail` 的字典，另有稀疏文件 `<前缀>.grm.sp` 的 `n_pairs`。
- **示例代码：**

```python
import tempfile
from pathlib import Path

out = Path(tempfile.mkdtemp())
from dnakit.gwas.gcta_p08_cli_functions import gcta_make_grm, gcta_make_sparse_grm

grm = gcta_make_grm(str(out / "g"), bfile="toy", verbose=False)["grm_prefix"]
sparse = gcta_make_sparse_grm(grm, str(out / "sp"), cutoff=0.05, verbose=False)
print(sparse["status"], sparse["n_pairs"], Path(sparse["outputs"]["grm_sp"]).name)
dense = gcta_make_sparse_grm(grm, str(out / "bk"), output_format="dense", verbose=False)
print(sorted(Path(p).name for k, p in dense["outputs"].items() if k != "log"))
```

- **示例输出：**

```text
completed 255 sp.grm.sp
['bk.grm.N.bin', 'bk.grm.bin', 'bk.grm.id']
```

## 4) `dnakit.gwas.gcta_inbreeding_coefficient()`

- **作用：** 估计每个个体的近交系数。
- **处理方法：** 运行 `gcta64 --autosome --ibc`，给出 Fhat1（方差标准化）、Fhat2（纯合过剩）、Fhat3（等位基因共享）。
- **API：** `dnakit.gwas.gcta_inbreeding_coefficient(bfile, out_prefix, *, keep=None, extract=None, maf=None, max_maf=None, threads=None, gcta_path=None, ..., verbose=True)`。
- **输入：** PLINK 前缀和新的输出前缀。
- **输出：** 包含 `status`、`backend`、`version`、`command`、`return_code`、`outputs`、`log_tail` 的字典，另有`summary`，含个体数和各估计的均值。
- **示例代码：**

```python
import tempfile
from pathlib import Path

out = Path(tempfile.mkdtemp())
from dnakit.gwas.gcta_p08_cli_functions import gcta_inbreeding_coefficient

result = gcta_inbreeding_coefficient("toy", str(out / "ibc"), verbose=False)
print(result["summary"]["individuals"], {k: round(v, 3) for k, v in result["summary"]["mean"].items()})
```

- **示例输出：**

```text
200 {'Fhat1': -0.002, 'Fhat2': -0.002, 'Fhat3': -0.002}
```

## 5) `dnakit.gwas.gcta_unify_grm()`

- **作用：** 把多个 GRM 限制到共同的个体并保持顺序一致。
- **处理方法：** 运行 `--mgrm list --unify-grm`；因 GCTA 自己决定输出位置和名称，函数在临时目录中运行并把结果移到 `<前缀>.<序号>.grm.*`。
- **API：** `dnakit.gwas.gcta_unify_grm(grm_prefixes, out_prefix, *, gcta_path=None, timeout_seconds=3600.0, verbose=True)`。
- **输入：** 至少两个 GRM 前缀和新的输出前缀。
- **输出：** 包含 `status`、`backend`、`version`、`command`、`return_code`、`outputs`、`log_tail` 的字典，另有`unified`：新的 GRM 前缀，以及 `common_individuals`。
- **示例代码：**

```python
import tempfile
from pathlib import Path

out = Path(tempfile.mkdtemp())
from dnakit.gwas.gcta_p08_cli_functions import gcta_make_grm, gcta_unify_grm

grm = gcta_make_grm(str(out / "g"), bfile="toy", verbose=False)["grm_prefix"]
subset = gcta_make_grm(str(out / "s"), bfile="toy", keep="keep150.txt", verbose=False)["grm_prefix"]
result = gcta_unify_grm([grm, subset], str(out / "un"), verbose=False)
print(result["status"], result["common_individuals"], [Path(p).name for p in result["unified"]])
```

- **示例输出：**

```text
completed 150 ['un.1', 'un.2']
```

## 6) `dnakit.gwas.gcta_reml()`

- **作用：** 用 GREML 估计方差组分和基于 SNP 的遗传力，可选 GxE 与遗传值预测。
- **处理方法：** 运行 `--reml`；`gxe_file` 加 `--gxe`，`predict='blup'` 加 `--reml-pred-rand`，`predict='cvblup'` 加 `--cvblup`。
- **API：** `dnakit.gwas.gcta_reml(phenotype_file, out_prefix, *, grm_prefix=None, grm_prefixes=None, phenotype_column=None, covariate_file=None, quantitative_covariate_file=None, gxe_file=None, ..., verbose=True)`。
- **输入：** 一个 GRM（`grm_prefix`）或多个 GRM（`grm_prefixes`）、表型文件和新的输出前缀。
- **输出：** 包含 `status`、`backend`、`version`、`command`、`return_code`、`outputs`、`log_tail` 的字典，另有`summary`（`variance_components` 与 `statistics`）；BLUP 时还有 `n_predictions`。
- **示例代码：**

```python
import tempfile
from pathlib import Path

out = Path(tempfile.mkdtemp())
from dnakit.gwas.gcta_p08_cli_functions import gcta_make_grm, gcta_reml

grm = gcta_make_grm(str(out / "g"), bfile="toy", verbose=False)["grm_prefix"]
result = gcta_reml("toy.phen", str(out / "r"), grm_prefix=grm, predict="blup", verbose=False)
for row in result["summary"]["variance_components"]:
    print(row["source"], round(row["variance"], 3), round(row["se"], 3))
print(result["summary"]["statistics"]["n"], result["n_predictions"])
gxe = gcta_reml("toy.phen", str(out / "gxe"), grm_prefix=grm, gxe_file="toy.env_class", verbose=False)
print([row["source"] for row in gxe["summary"]["variance_components"]][:3])
```

- **示例输出：**

```text
V(G) 1.877 2.179
V(e) 4.319 2.169
Vp 6.196 0.623
V(G)/Vp 0.303 0.348
200 200
['V(G)', 'V(GxE)', 'V(e)']
```

## 7) `dnakit.gwas.gcta_blup_snp()`

- **作用：** 把个体遗传值的 BLUP 转换成每个 SNP 的 BLUP 效应。
- **处理方法：** 运行 `--blup-snp`，输入是 `gcta_reml(predict='blup')` 写出的 `.indi.blp`。
- **API：** `dnakit.gwas.gcta_blup_snp(bfile, indi_blp_file, out_prefix, *, autosome=False, keep=None, extract=None, threads=None, gcta_path=None, timeout_seconds=3600.0, verbose=True)`。
- **输入：** PLINK 前缀、`.indi.blp` 文件和新的输出前缀（含性染色体时设 `autosome=True`）。
- **输出：** 包含 `status`、`backend`、`version`、`command`、`return_code`、`outputs`、`log_tail` 的字典，另有`n_snps`；效应在 `<前缀>.snp.blp`。
- **示例代码：**

```python
import tempfile
from pathlib import Path

out = Path(tempfile.mkdtemp())
from dnakit.gwas.gcta_p08_cli_functions import gcta_blup_snp, gcta_make_grm, gcta_reml

grm = gcta_make_grm(str(out / "g"), bfile="toy", verbose=False)["grm_prefix"]
blup = gcta_reml("toy.phen", str(out / "r"), grm_prefix=grm, predict="blup", verbose=False)
result = gcta_blup_snp("toy", blup["outputs"]["indi_blp"], str(out / "bs"), autosome=True, verbose=False)
print(result["status"], result["n_snps"])
print(gcta_blup_snp("toy", blup["outputs"]["indi_blp"], str(out / "bs2"), verbose=False)["status"])
```

- **示例输出：**

```text
completed 2500
failed
```

## 8) `dnakit.gwas.gcta_reml_bivar()`

- **作用：** 双变量 GREML：估计两个性状之间的遗传相关 rG。
- **处理方法：** 运行 `--reml-bivar`。
- **API：** `dnakit.gwas.gcta_reml_bivar(phenotype_file, out_prefix, *, grm_prefix=None, grm_prefixes=None, trait_columns=None, covariate_file=None, quantitative_covariate_file=None, no_constrain=False, ..., verbose=True)`。
- **输入：** GRM、含两个性状的表型文件和新的输出前缀。
- **输出：** 包含 `status`、`backend`、`version`、`command`、`return_code`、`outputs`、`log_tail` 的字典，另有`summary`，含包括 `rG` 在内的方差组分。
- **示例代码：**

```python
import tempfile
from pathlib import Path

out = Path(tempfile.mkdtemp())
from dnakit.gwas.gcta_p08_cli_functions import gcta_make_grm, gcta_reml_bivar

grm = gcta_make_grm(str(out / "g"), bfile="toy", verbose=False)["grm_prefix"]
result = gcta_reml_bivar("toy2.phen", str(out / "bv"), grm_prefix=grm, verbose=False)
print([row["source"] for row in result["summary"]["variance_components"]])
```

- **示例输出：**

```text
['V(G)_tr1', 'V(G)_tr2', 'C(G)_tr12', 'V(e)_tr1', 'V(e)_tr2', 'C(e)_tr12', 'Vp_tr1', 'Vp_tr2', 'V(G)/Vp_tr1', 'V(G)/Vp_tr2', 'rG']
```

## 9) `dnakit.gwas.gcta_he_regression()`

- **作用：** Haseman-Elston 回归：估计遗传力或（双变量时）遗传相关。
- **处理方法：** 运行 `--HEreg` 或 `--HEreg-bivar`，给出 HE-CP（交叉积）和 HE-SD（平方差）两张表。
- **API：** `dnakit.gwas.gcta_he_regression(phenotype_file, out_prefix, *, grm_prefix=None, grm_prefixes=None, bivariate=False, trait_columns=(1, 2), phenotype_column=None, keep=None, ..., verbose=True)`。
- **输入：** GRM、表型文件和新的输出前缀；两个性状时设 `bivariate=True`。
- **输出：** 包含 `status`、`backend`、`version`、`command`、`return_code`、`outputs`、`log_tail` 的字典，另有`summary`，含 `HE-CP` 与 `HE-SD` 两张表。
- **示例代码：**

```python
import tempfile
from pathlib import Path

out = Path(tempfile.mkdtemp())
from dnakit.gwas.gcta_p08_cli_functions import gcta_he_regression, gcta_make_grm

grm = gcta_make_grm(str(out / "g"), bfile="toy", verbose=False)["grm_prefix"]
he = gcta_he_regression("toy.phen", str(out / "he"), grm_prefix=grm, verbose=False)
print(sorted(he["summary"]), [row["coefficient"] for row in he["summary"]["HE-CP"]])
pair = gcta_he_regression("toy2.phen", str(out / "hb"), grm_prefix=grm, bivariate=True, verbose=False)
print([row["coefficient"] for row in pair["summary"]["HE-CP"]][-3:])
```

- **示例输出：**

```text
['HE-CP', 'HE-SD'] ['Intercept', 'V(G)/Vp']
['rG', 'N_tr1', 'N_tr2']
```

## 10) `dnakit.gwas.gcta_fastgwa()`

- **作用：** fastGWA 关联检验（线性混合模型、精确模型、二分类 GLMM、线性回归及 GxE）。
- **处理方法：** `model` 取 `mlm`、`mlm_exact`、`mlm_binary` 或 `lr`；`envir_file` 给出时做 fastGWA-GE；`model_only=True` 只保存模型。
- **API：** `dnakit.gwas.gcta_fastgwa(phenotype_file, out_prefix, *, bfile=None, mbfile=None, model='mlm', sparse_grm_prefix=None, covariate_file=None, quantitative_covariate_file=None, ..., verbose=True)`。
- **输入：** 表型文件、输出前缀、基因型，混合模型还需稀疏 GRM。
- **输出：** 包含 `status`、`backend`、`version`、`command`、`return_code`、`outputs`、`log_tail` 的字典，另有`<前缀>.fastGWA` 的 `summary`（行数、列名和 p 最小的行）。
- **示例代码：**

```python
import tempfile
from pathlib import Path

out = Path(tempfile.mkdtemp())
from dnakit.gwas.gcta_p08_cli_functions import gcta_fastgwa, gcta_make_grm, gcta_make_sparse_grm

grm = gcta_make_grm(str(out / "g"), bfile="toy", verbose=False)["grm_prefix"]
sparse = gcta_make_sparse_grm(grm, str(out / "sp"), verbose=False)["grm_prefix"]
linear = gcta_fastgwa("toy.phen", str(out / "lr"), bfile="toy", model="lr", verbose=False)
print(linear["summary"]["rows"], linear["summary"]["top_rows"][0]["SNP"])
mixed = gcta_fastgwa("toy.phen", str(out / "mlm"), bfile="toy", sparse_grm_prefix=sparse, verbose=False)
print(mixed["summary"]["columns"])
ge = gcta_fastgwa("toy.phen", str(out / "ge"), bfile="toy", model="lr", envir_file="toy.envir", verbose=False)
print(ge["summary"]["columns"][7:10])
binary = gcta_fastgwa("toy.bin.phen", str(out / "bin"), bfile="toy", model="mlm_binary", sparse_grm_prefix=sparse, quantitative_covariate_file="toy.qcovar", verbose=False)
print(binary["summary"]["columns"][-1], binary["summary"]["rows"])
```

- **示例输出：**

```text
2600 rs223
['CHR', 'SNP', 'POS', 'A1', 'A2', 'N', 'AF1', 'BETA', 'SE', 'P']
['BETA_G', 'BETA_G_by_E', 'SE_G']
CONVERGE 2600
```

## 11) `dnakit.gwas.gcta_mlma()`

- **作用：** 混合线性模型关联分析，可选留一条染色体（LOCO）。
- **处理方法：** 运行 `--mlma`、`--mlma-loco` 或 `--mlma-subtract-grm`。
- **API：** `dnakit.gwas.gcta_mlma(bfile, phenotype_file, out_prefix, *, grm_prefix=None, grm_prefixes=None, loco=False, subtract_grm_prefix=None, covariate_file=None, quantitative_covariate_file=None, ..., verbose=True)`。
- **输入：** 基因型、表型文件、新的输出前缀，以及 GRM（或 `loco=True`）。
- **输出：** 包含 `status`、`backend`、`version`、`command`、`return_code`、`outputs`、`log_tail` 的字典，另有`<前缀>.mlma` 或 `.loco.mlma` 的 `summary`。
- **示例代码：**

```python
import tempfile
from pathlib import Path

out = Path(tempfile.mkdtemp())
from dnakit.gwas.gcta_p08_cli_functions import gcta_make_grm, gcta_mlma

grm = gcta_make_grm(str(out / "g"), bfile="toy", verbose=False)["grm_prefix"]
plain = gcta_mlma("toy", "toy.phen", str(out / "m"), grm_prefix=grm, verbose=False)
print(plain["summary"]["rows"], plain["summary"]["columns"], plain["summary"]["top_rows"][0]["SNP"])
loco = gcta_mlma("toy", "toy.phen", str(out / "l"), loco=True, verbose=False)
print(Path(loco["outputs"]["mlma"]).name)
```

- **示例输出：**

```text
2600 ['Chr', 'SNP', 'bp', 'A1', 'A2', 'Freq', 'b', 'se', 'p'] rs223
l.loco.mlma
```

## 12) `dnakit.gwas.gcta_fastgwa_burden()`

- **作用：** fastGWA-BB：二分类性状的基因级稀有变异负荷检验。
- **处理方法：** 运行 `--burden`，使用 `gcta_fastgwa(model='mlm_binary', model_only=True)` 保存的模型。
- **API：** `dnakit.gwas.gcta_fastgwa_burden(bfile, model_prefix, set_list_file, out_prefix, *, max_maf=None, maf=None, extract=None, geno=None, threads=None, gcta_path=None, ..., verbose=True)`。
- **输入：** 基因型、已保存的模型前缀、基因列表和新的输出前缀。
- **输出：** 包含 `status`、`backend`、`version`、`command`、`return_code`、`outputs`、`log_tail` 的字典，另有`<前缀>.fastGWA` 的 `summary`，每个基因一行。
- **示例代码：**

```python
import tempfile
from pathlib import Path

out = Path(tempfile.mkdtemp())
from dnakit.gwas.gcta_p08_cli_functions import gcta_fastgwa, gcta_fastgwa_burden, gcta_make_grm, gcta_make_sparse_grm

grm = gcta_make_grm(str(out / "g"), bfile="toy", verbose=False)["grm_prefix"]
sparse = gcta_make_sparse_grm(grm, str(out / "sp"), verbose=False)["grm_prefix"]
model = gcta_fastgwa("toy.bin.phen", str(out / "mdl"), bfile="toy", model="mlm_binary", sparse_grm_prefix=sparse, quantitative_covariate_file="toy.qcovar", model_only=True, verbose=False)
result = gcta_fastgwa_burden("toy", model["model_prefix"], "genes_burden.txt", str(out / "bur"), max_maf=0.5, maf=0.0001, verbose=False)
print(result["summary"]["rows"], [row["GENE"] for row in result["summary"]["top_rows"]], [row["VAR_N"] for row in sorted(result["summary"]["top_rows"], key=lambda r: r["GENE"])])
```

- **示例输出：**

```text
3 ['G2', 'G3', 'G1'] [100, 201, 300]
```

## 13) `dnakit.gwas.gcta_acat()`

- **作用：** 由单变异 p 值做基因级稀有变异 ACAT-V 检验。
- **处理方法：** 运行 `--acat`，输入是 fastGWA 关联表。
- **API：** `dnakit.gwas.gcta_acat(snp_list_file, gene_list_file, output_path, *, max_maf=None, min_mac=None, window_kb=None, gcta_path=None, timeout_seconds=3600.0, verbose=True)`。
- **输入：** fastGWA 关联表、基因列表（染色体、起点、终点、基因名）和新的输出路径。
- **输出：** 包含 `status`、`backend`、`version`、`command`、`return_code`、`outputs`、`log_tail` 的字典，另有结果表 `summary`（CHR、GENE、START、END、SNP_NUM、P_ACAT）。
- **示例代码：**

```python
import tempfile
from pathlib import Path

out = Path(tempfile.mkdtemp())
from dnakit.gwas.gcta_p08_cli_functions import gcta_acat, gcta_fastgwa, gcta_make_grm, gcta_make_sparse_grm

grm = gcta_make_grm(str(out / "g"), bfile="toy", verbose=False)["grm_prefix"]
sparse = gcta_make_sparse_grm(grm, str(out / "sp"), verbose=False)["grm_prefix"]
single = gcta_fastgwa("toy.bin.phen", str(out / "single"), bfile="toy", model="mlm_binary", sparse_grm_prefix=sparse, quantitative_covariate_file="toy.qcovar", verbose=False)
result = gcta_acat(single["outputs"]["fastgwa"], "genes_acat.txt", str(out / "acat.res"), max_maf=0.5, min_mac=5, verbose=False)
print(sorted((row["GENE"], row["SNP_NUM"]) for row in result["summary"]["top_rows"]))
```

- **示例输出：**

```text
[('G1', 100), ('G2', 201), ('G3', 300)]
```

## 14) `dnakit.gwas.gcta_cojo()`

- **作用：** 基于 GWAS 汇总统计的条件与联合多 SNP 分析（COJO）。
- **处理方法：** `mode='slct'` 逐步选择独立信号（`top_snps` 选固定个数），`'joint'` 联合拟合，`'cond'` 以指定 SNP 为条件。
- **API：** `dnakit.gwas.gcta_cojo(bfile, summary_file, out_prefix, *, mode='slct', snp_list_file=None, top_snps=None, p_threshold=None, window_kb=None, collinearity=None, ..., verbose=True)`。
- **输入：** 基因型（LD 参考）、COJO `.ma` 汇总统计文件、新的输出前缀和 `mode`。
- **输出：** 包含 `status`、`backend`、`version`、`command`、`return_code`、`outputs`、`log_tail` 的字典，另有`.jma`（或 `.cma`）表的 `summary`。
- **示例代码：**

```python
import tempfile
from pathlib import Path

out = Path(tempfile.mkdtemp())
from dnakit.gwas.gcta_p08_cli_functions import gcta_cojo, gcta_make_grm, gcta_mlma

grm = gcta_make_grm(str(out / "g"), bfile="toy", verbose=False)["grm_prefix"]
ml = gcta_mlma("toy", "toy.phen", str(out / "ml"), grm_prefix=grm, verbose=False)
rows = [line.split() for line in Path(ml["outputs"]["mlma"]).read_text().splitlines()[1:]]
(out / "ml.ma").write_text("SNP A1 A2 freq b se p N\n" + "\n".join(" ".join([r[1], r[3], r[4], r[5], r[6], r[7], r[8], "200"]) for r in rows) + "\n")
ma = str(out / "ml.ma")
selected = gcta_cojo("toy", ma, str(out / "s"), chromosome=1, p_threshold=1e-2, verbose=False)
print(selected["summary"]["rows"], sorted(selected["outputs"])[:3])
top = gcta_cojo("toy", ma, str(out / "t"), chromosome=1, top_snps=3, verbose=False)
print(top["summary"]["rows"])
cond = gcta_cojo("toy", ma, str(out / "c"), mode="cond", snp_list_file="cond.snplist", chromosome=1, verbose=False)
print(cond["summary"]["rows"], "given" in cond["outputs"])
```

- **示例输出：**

```text
7 ['cma', 'jma', 'ldr']
3
1299 True
```

## 15) `dnakit.gwas.gcta_sblup()`

- **作用：** 由汇总统计计算 SBLUP 风险预测的 SNP 效应。
- **处理方法：** 运行 `--cojo-sblup`；`shrinkage` 为 λ = m(1/h²−1)。
- **API：** `dnakit.gwas.gcta_sblup(bfile, summary_file, out_prefix, shrinkage, *, window_kb=200, keep=None, extract=None, chromosome=None, maf=None, threads=None, ..., verbose=True)`。
- **输入：** 基因型、`.ma` 汇总统计文件、新的输出前缀和 `shrinkage`。
- **输出：** 包含 `status`、`backend`、`version`、`command`、`return_code`、`outputs`、`log_tail` 的字典，另有`n_snps`；效应在 `<前缀>.sblup.cojo`。
- **示例代码：**

```python
import tempfile
from pathlib import Path

out = Path(tempfile.mkdtemp())
from dnakit.gwas.gcta_p08_cli_functions import gcta_make_grm, gcta_mlma, gcta_sblup

grm = gcta_make_grm(str(out / "g"), bfile="toy", verbose=False)["grm_prefix"]
ml = gcta_mlma("toy", "toy.phen", str(out / "ml"), grm_prefix=grm, verbose=False)
rows = [line.split() for line in Path(ml["outputs"]["mlma"]).read_text().splitlines()[1:]]
(out / "ml.ma").write_text("SNP A1 A2 freq b se p N\n" + "\n".join(" ".join([r[1], r[3], r[4], r[5], r[6], r[7], r[8], "200"]) for r in rows) + "\n")
result = gcta_sblup("toy", str(out / "ml.ma"), str(out / "sb"), 2000, chromosome=1, verbose=False)
print(result["status"], result["n_snps"])
```

- **示例输出：**

```text
completed 1300
```

## 16) `dnakit.gwas.gcta_mtcojo()`

- **作用：** 多性状条件与联合分析：以其它性状为条件做目标性状的 GWAS。
- **处理方法：** 运行 `--mtcojo-file`；协变量性状对目标性状的效应 bxy 由 GSMR 估计，或来自 `bxy_file`。
- **API：** `dnakit.gwas.gcta_mtcojo(summary_list_file, out_prefix, ref_ld_dir, w_ld_dir, *, bfile=None, mbfile=None, bxy_file=None, diff_freq=None, gwas_thresh=None, clump_r2=None, ..., verbose=True)`。
- **输入：** 性状及汇总统计文件列表、LDSC 格式 LD score 目录、基因型（LD 参考）和新的输出前缀。
- **输出：** 包含 `status`、`backend`、`version`、`command`、`return_code`、`outputs`、`log_tail` 的字典，另有`<前缀>.mtcojo.cma` 的 `summary`。
- **示例代码：**

```python
import tempfile
from pathlib import Path

out = Path(tempfile.mkdtemp())
from dnakit.gwas.gcta_p08_cli_functions import gcta_make_grm, gcta_mlma, gcta_mtcojo
import gzip

grm = gcta_make_grm(str(out / "g"), bfile="toy", verbose=False)["grm_prefix"]
ml = gcta_mlma("toy", "toy.phen", str(out / "ml"), grm_prefix=grm, verbose=False)
rows = [line.split() for line in Path(ml["outputs"]["mlma"]).read_text().splitlines()[1:]]
(out / "ml.ma").write_text("SNP A1 A2 freq b se p N\n" + "\n".join(" ".join([r[1], r[3], r[4], r[5], r[6], r[7], r[8], "200"]) for r in rows) + "\n")
lines = (out / "ml.ma").read_text().splitlines()
(out / "exp.ma").write_text(lines[0] + "\n" + "\n".join(" ".join(l.split()[:4] + [str(float(l.split()[4]) * 0.5 + 0.01)] + l.split()[5:]) for l in lines[1:]) + "\n")
(out / "mt.list").write_text(f"t2d {out / 'ml.ma'}\nbmi {out / 'exp.ma'}\n")
(out / "bxy.txt").write_text("bmi 0.5\n")
ld = out / "ldsc"
ld.mkdir()
bim = [line.split() for line in Path("toy.bim").read_text().splitlines()]
for chromosome in range(1, 23):
    rows = [row for row in bim if row[0] == str(chromosome)]
    with gzip.open(ld / f"{chromosome}.l2.ldscore.gz", "wt") as handle:
        handle.write("CHR\tSNP\tBP\tL2\n" + "".join(f"{r[0]}\t{r[1]}\t{r[3]}\t1.5\n" for r in rows))
    (ld / f"{chromosome}.l2.M_5_50").write_text(f"{len(rows)}\n")
result = gcta_mtcojo(str(out / "mt.list"), str(out / "mt"), str(ld), str(ld), bfile="toy", bxy_file=str(out / "bxy.txt"), verbose=False)
print(result["status"], result["summary"]["rows"], result["summary"]["columns"][-3:])
```

- **示例输出：**

```text
completed 2600 ['bC', 'bC_se', 'bC_pval']
```

## 17) `dnakit.gwas.gcta_fastbat()`

- **作用：** 基于 GWAS 汇总统计和 LD 参考的集合/基因/区段关联检验（fastBAT）。
- **处理方法：** 运行 `--fastBAT`，集合定义三选一：基因列表、SNP 集合列表或区段长度。
- **API：** `dnakit.gwas.gcta_fastbat(bfile, summary_file, out_prefix, *, gene_list_file=None, set_list_file=None, segment_kb=None, gene_window_kb=None, ld_cutoff=None, diff_freq=None, ..., verbose=True)`。
- **输入：** 基因型、`.ma` 汇总统计文件、新的输出前缀，以及 `gene_list_file`、`set_list_file`、`segment_kb` 之一。
- **输出：** 包含 `status`、`backend`、`version`、`command`、`return_code`、`outputs`、`log_tail` 的字典，另有`.gene.fastbat`、`.fastbat` 或 `.seg.fastbat` 的 `summary`。
- **示例代码：**

```python
import tempfile
from pathlib import Path

out = Path(tempfile.mkdtemp())
from dnakit.gwas.gcta_p08_cli_functions import gcta_fastbat, gcta_make_grm, gcta_mlma

grm = gcta_make_grm(str(out / "g"), bfile="toy", verbose=False)["grm_prefix"]
ml = gcta_mlma("toy", "toy.phen", str(out / "ml"), grm_prefix=grm, verbose=False)
rows = [line.split() for line in Path(ml["outputs"]["mlma"]).read_text().splitlines()[1:]]
(out / "ml.ma").write_text("SNP A1 A2 freq b se p N\n" + "\n".join(" ".join([r[1], r[3], r[4], r[5], r[6], r[7], r[8], "200"]) for r in rows) + "\n")
ma = str(out / "ml.ma")
gene = gcta_fastbat("toy", ma, str(out / "gene"), gene_list_file="genes.txt", verbose=False)
print(Path(gene["outputs"]["fastbat"]).name, gene["summary"]["rows"], gene["summary"]["columns"][:5])
segment = gcta_fastbat("toy", ma, str(out / "seg"), segment_kb=100, verbose=False)
print(Path(segment["outputs"]["fastbat"]).name, segment["summary"]["rows"])
sets = gcta_fastbat("toy", ma, str(out / "sets"), set_list_file="sets.txt", verbose=False)
print(Path(sets["outputs"]["fastbat"]).name, sets["summary"]["rows"])
```

- **示例输出：**

```text
gene.gene.fastbat 3 ['Gene', 'Chr', 'Start', 'End', 'No.SNPs']
seg.seg.fastbat 26
sets.fastbat 2
```

## 18) `dnakit.gwas.gcta_write_gene_snpsets()`

- **作用：** 构建并写出基因级检验所用的每个基因的 SNP 集合。
- **处理方法：** 运行 `--fastBAT-write-snpset`（`method='fastbat'`）或 `--mBAT-write-snpset`（`method='mbat'`）；基因检验本身也会运行。
- **API：** `dnakit.gwas.gcta_write_gene_snpsets(bfile, summary_file, gene_list_file, out_prefix, *, method='fastbat', gene_window_kb=None, diff_freq=None, keep=None, extract=None, chromosome=None, ..., verbose=True)`。
- **输入：** 基因型、`.ma` 汇总统计文件、基因列表和新的输出前缀。
- **输出：** 包含 `status`、`backend`、`version`、`command`、`return_code`、`outputs`、`log_tail` 的字典，另有`snpset`（集合文件路径）与 `genes`。
- **示例代码：**

```python
import tempfile
from pathlib import Path

out = Path(tempfile.mkdtemp())
from dnakit.gwas.gcta_p08_cli_functions import gcta_fastbat, gcta_make_grm, gcta_mlma, gcta_write_gene_snpsets

grm = gcta_make_grm(str(out / "g"), bfile="toy", verbose=False)["grm_prefix"]
ml = gcta_mlma("toy", "toy.phen", str(out / "ml"), grm_prefix=grm, verbose=False)
rows = [line.split() for line in Path(ml["outputs"]["mlma"]).read_text().splitlines()[1:]]
(out / "ml.ma").write_text("SNP A1 A2 freq b se p N\n" + "\n".join(" ".join([r[1], r[3], r[4], r[5], r[6], r[7], r[8], "200"]) for r in rows) + "\n")
ma = str(out / "ml.ma")
made = gcta_write_gene_snpsets("toy", ma, "genes.txt", str(out / "ss"), verbose=False)
print(Path(made["snpset"]).name, made["genes"])
reuse = gcta_fastbat("toy", ma, str(out / "again"), set_list_file=made["snpset"], verbose=False)
print(reuse["summary"]["rows"])
mbat = gcta_write_gene_snpsets("toy", ma, "genes.txt", str(out / "ms"), method="mbat", verbose=False)
print(Path(mbat["snpset"]).name, mbat["genes"])
```

- **示例输出：**

```text
ss.gene.snpset 3
3
ms.gene.snpset.mbat 3
```

## 19) `dnakit.gwas.gcta_mbat_combo()`

- **作用：** mBAT-combo 基因级检验：存在效应相互掩盖时仍有较高功效。
- **处理方法：** 运行 `--mBAT-combo`，结合 mBAT（LD 矩阵主成分）与 fastBAT 两种检验。
- **API：** `dnakit.gwas.gcta_mbat_combo(bfile, summary_file, gene_list_file, out_prefix, *, gene_window_kb=None, svd_gamma=None, ld_cutoff=None, print_all_p=False, diff_freq=None, keep=None, ..., verbose=True)`。
- **输入：** 基因型、`.ma` 汇总统计文件、基因列表和新的输出前缀。
- **输出：** 包含 `status`、`backend`、`version`、`command`、`return_code`、`outputs`、`log_tail` 的字典，另有`<前缀>.gene.assoc.mbat` 的 `summary`。
- **示例代码：**

```python
import tempfile
from pathlib import Path

out = Path(tempfile.mkdtemp())
from dnakit.gwas.gcta_p08_cli_functions import gcta_make_grm, gcta_mbat_combo, gcta_mlma

grm = gcta_make_grm(str(out / "g"), bfile="toy", verbose=False)["grm_prefix"]
ml = gcta_mlma("toy", "toy.phen", str(out / "ml"), grm_prefix=grm, verbose=False)
rows = [line.split() for line in Path(ml["outputs"]["mlma"]).read_text().splitlines()[1:]]
(out / "ml.ma").write_text("SNP A1 A2 freq b se p N\n" + "\n".join(" ".join([r[1], r[3], r[4], r[5], r[6], r[7], r[8], "200"]) for r in rows) + "\n")
result = gcta_mbat_combo("toy", str(out / "ml.ma"), "genes.txt", str(out / "mb"), print_all_p=True, verbose=False)
print(result["summary"]["rows"], [c for c in result["summary"]["columns"] if c.startswith("P_")])
```

- **示例输出：**

```text
3 ['P_mBATcombo', 'P_mBAT', 'P_fastBAT']
```

## 20) `dnakit.gwas.gcta_gsmr()`

- **作用：** 广义汇总数据孟德尔随机化（GSMR）：暴露对结局的因果效应。
- **处理方法：** 运行 `--gsmr-file`，对暴露的 SNP 工具变量做 clumping 和 HEIDI 离群值过滤。
- **API：** `dnakit.gwas.gcta_gsmr(exposure_list_file, outcome_list_file, out_prefix, *, bfile=None, mbfile=None, direction=0, gsmr2=False, diff_freq=None, gwas_thresh=None, ..., verbose=True)`。
- **输入：** 基因型（LD 参考）、暴露与结局列表文件和新的输出前缀。
- **输出：** 包含 `status`、`backend`、`version`、`command`、`return_code`、`outputs`、`log_tail` 的字典，另有`<前缀>.gsmr` 的 `summary`（exposure、outcome、bxy、se、p、nsnp）。
- **示例代码：**

```python
import tempfile
from pathlib import Path

out = Path(tempfile.mkdtemp())
from dnakit.gwas.gcta_p08_cli_functions import gcta_gsmr, gcta_make_grm, gcta_mlma

grm = gcta_make_grm(str(out / "g"), bfile="toy", verbose=False)["grm_prefix"]
ml = gcta_mlma("toy", "toy.phen", str(out / "ml"), grm_prefix=grm, verbose=False)
rows = [line.split() for line in Path(ml["outputs"]["mlma"]).read_text().splitlines()[1:]]
(out / "ml.ma").write_text("SNP A1 A2 freq b se p N\n" + "\n".join(" ".join([r[1], r[3], r[4], r[5], r[6], r[7], r[8], "200"]) for r in rows) + "\n")
lines = (out / "ml.ma").read_text().splitlines()
(out / "exp.ma").write_text(lines[0] + "\n" + "\n".join(" ".join(l.split()[:4] + [str(float(l.split()[4]) * 0.5 + 0.01)] + l.split()[5:]) for l in lines[1:]) + "\n")
(out / "expo.txt").write_text(f"bmi {out / 'exp.ma'}\n")
(out / "outc.txt").write_text(f"t2d {out / 'ml.ma'}\n")
result = gcta_gsmr(str(out / "expo.txt"), str(out / "outc.txt"), str(out / "gs"), bfile="toy", snp_min=3, gwas_thresh=0.05, clump_r2=0.5, heidi_thresh=0.0, verbose=False)
row = result["summary"]["top_rows"][0]
print(row["Exposure"], row["Outcome"], row["nsnp"], row["bxy"] > 0)
```

- **示例输出：**

```text
bmi t2d 3 True
```

## 21) `dnakit.gwas.gcta_ld_friends()`

- **作用：** 为每个目标 SNP 用简单回归寻找与其显著连锁的 SNP。
- **处理方法：** 运行 `--ld`，给出每个目标的 LD 区块范围、SNP 数和 r² 统计。
- **API：** `dnakit.gwas.gcta_ld_friends(bfile, snp_list_file, out_prefix, *, window_kb=None, significance=None, keep=None, extract=None, threads=None, gcta_path=None, ..., verbose=True)`。
- **输入：** 基因型、目标 SNP 列表和新的输出前缀。
- **输出：** 包含 `status`、`backend`、`version`、`command`、`return_code`、`outputs`、`log_tail` 的字典，另有`<前缀>.rsq.ld` 的 `summary`；另有 `.r.ld` 与 `.snp.ld`。
- **示例代码：**

```python
import tempfile
from pathlib import Path

out = Path(tempfile.mkdtemp())
from dnakit.gwas.gcta_p08_cli_functions import gcta_ld_friends

result = gcta_ld_friends("toy", "ld.snplist", str(out / "ldf"), window_kb=5000, significance=0.05, verbose=False)
print([row["target_SNP"] for row in result["summary"]["first_rows"]], sorted(k for k in result["outputs"] if k != "log"))
```

- **示例输出：**

```text
['rs5', 'rs100'] ['r', 'rsq', 'snp']
```

## 22) `dnakit.gwas.gcta_segment_ld_score()`

- **作用：** 计算基于区段的 LD score，用于 GREML-LDMS 的 SNP 分层。
- **处理方法：** 运行 `--autosome --ld-score-region`；区段默认 200 kb，相邻区段重叠一半。
- **API：** `dnakit.gwas.gcta_segment_ld_score(bfile, out_prefix, *, segment_kb=200, window_kb=None, rsq_cutoff=None, adjust=False, keep=None, extract=None, ..., verbose=True)`。
- **输入：** 基因型和新的输出前缀；可选 `segment_kb`。
- **输出：** 包含 `status`、`backend`、`version`、`command`、`return_code`、`outputs`、`log_tail` 的字典，另有`<前缀>.score.ld` 的 `summary`。
- **示例代码：**

```python
import tempfile
from pathlib import Path

out = Path(tempfile.mkdtemp())
from dnakit.gwas.gcta_p08_cli_functions import gcta_segment_ld_score

result = gcta_segment_ld_score("toy", str(out / "seg"), verbose=False)
print(result["summary"]["rows"], result["summary"]["columns"])
```

- **示例输出：**

```text
2500 ['SNP', 'chr', 'bp', 'freq', 'mean_rsq', 'snp_num', 'max_rsq', 'ldscore_SNP', 'ldscore_region']
```

## 23) `dnakit.gwas.gcta_multi_set_ld_score()`

- **作用：** 计算每个 SNP 相对多个 SNP 集合的 LD score，用于多组分 LD score 回归。
- **处理方法：** 运行 `--autosome --ld-score-multi`，函数把集合文件列表写入临时文件。
- **API：** `dnakit.gwas.gcta_multi_set_ld_score(bfile, snp_set_files, out_prefix, *, window_kb=None, rsq_cutoff=None, adjust=False, keep=None, extract=None, maf=None, ..., verbose=True)`。
- **输入：** 基因型、SNP 集合文件列表和新的输出前缀。
- **输出：** 包含 `status`、`backend`、`version`、`command`、`return_code`、`outputs`、`log_tail` 的字典，另有`n_sets` 与得分表的 `summary`。
- **示例代码：**

```python
import tempfile
from pathlib import Path

out = Path(tempfile.mkdtemp())
from dnakit.gwas.gcta_p08_cli_functions import gcta_multi_set_ld_score

result = gcta_multi_set_ld_score("toy", ["set1.snps", "set2.snps"], str(out / "multi"), window_kb=1000, verbose=False)
print(result["n_sets"], result["summary"]["rows"], Path(result["outputs"]["score_ld"]).name)
```

- **示例输出：**

```text
2 2500 multi.multi_score.ld
```

## 24) `dnakit.gwas.gcta_pc_loading()`

- **作用：** 由参考样本的主成分计算每个 SNP 的载荷。
- **处理方法：** 运行 `--pc-loading`；需要 GCTA `--pca` 写出的 `.eigenvec` 和 `.eigenval`。
- **API：** `dnakit.gwas.gcta_pc_loading(bfile, pca_prefix, out_prefix, *, autosome=True, keep=None, extract=None, maf=None, threads=None, gcta_path=None, ..., verbose=True)`。
- **输入：** 基因型、PCA 文件前缀和新的输出前缀。
- **输出：** 包含 `status`、`backend`、`version`、`command`、`return_code`、`outputs`、`log_tail` 的字典，另有`n_snps`；载荷在 `<前缀>.pcl`。
- **示例代码：**

```python
import tempfile
from pathlib import Path

out = Path(tempfile.mkdtemp())
from dnakit.gwas.gcta_p08_cli_functions import gcta_make_grm, gcta_pc_loading
import subprocess

grm = gcta_make_grm(str(out / "g"), bfile="toy", verbose=False)["grm_prefix"]
subprocess.run(["gcta64", "--grm", grm, "--pca", "5", "--out", str(out / "pc")], check=True, capture_output=True)
result = gcta_pc_loading("toy", str(out / "pc"), str(out / "pl"), verbose=False)
print(result["status"], result["n_snps"])
```

- **示例输出：**

```text
completed 2500
```

## 25) `dnakit.gwas.gcta_project_loading()`

- **作用：** 把样本投影到参考 SNP 载荷上得到主成分。
- **处理方法：** 运行 `--project-loading`，仅使用与参考共有的 SNP。
- **API：** `dnakit.gwas.gcta_project_loading(bfile, loading_prefix, out_prefix, *, n_components=20, autosome=True, update_freq_file=None, keep=None, extract=None, threads=None, ..., verbose=True)`。
- **输入：** 基因型、`.pcl` 文件前缀和新的输出前缀。
- **输出：** 包含 `status`、`backend`、`version`、`command`、`return_code`、`outputs`、`log_tail` 的字典，另有`n_individuals`；主成分在 `<前缀>.proj.eigenvec`。
- **示例代码：**

```python
import tempfile
from pathlib import Path

out = Path(tempfile.mkdtemp())
from dnakit.gwas.gcta_p08_cli_functions import gcta_make_grm, gcta_pc_loading, gcta_project_loading
import subprocess

grm = gcta_make_grm(str(out / "g"), bfile="toy", verbose=False)["grm_prefix"]
subprocess.run(["gcta64", "--grm", grm, "--pca", "5", "--out", str(out / "pc")], check=True, capture_output=True)
gcta_pc_loading("toy", str(out / "pc"), str(out / "pl"), verbose=False)
result = gcta_project_loading("toy", str(out / "pl"), str(out / "proj"), n_components=5, verbose=False)
reference = [float(line.split()[2]) for line in Path(str(out / "pc") + ".eigenvec").read_text().splitlines()]
projected = [float(line.split()[2]) for line in Path(result["outputs"]["eigenvec"]).read_text().splitlines()]
mx, my = sum(reference) / 200, sum(projected) / 200
cov = sum((a - mx) * (b - my) for a, b in zip(reference, projected))
corr = cov / (sum((a - mx) ** 2 for a in reference) ** 0.5 * sum((b - my) ** 2 for b in projected) ** 0.5)
print(result["n_individuals"], round(abs(corr), 3))
```

- **示例输出：**

```text
200 1.0
```

## 26) `dnakit.gwas.gcta_fst()`

- **作用：** 按亚群体计算每个 SNP 的 Fst（Weir 1996）。
- **处理方法：** 运行 `--fst --sub-popu`。
- **API：** `dnakit.gwas.gcta_fst(bfile, subpopulation_file, out_prefix, *, keep=None, extract=None, chromosome=None, maf=None, threads=None, gcta_path=None, ..., verbose=True)`。
- **输入：** 基因型、亚群体文件（FID、IID、标签）和新的输出前缀。
- **输出：** 包含 `status`、`backend`、`version`、`command`、`return_code`、`outputs`、`log_tail` 的字典，另有`summary`，含列名和 `mean_fst`。
- **示例代码：**

```python
import tempfile
from pathlib import Path

out = Path(tempfile.mkdtemp())
from dnakit.gwas.gcta_p08_cli_functions import gcta_fst

result = gcta_fst("toy", "toy.subpop", str(out / "fst"), verbose=False)
print(result["summary"]["rows"], result["summary"]["columns"][4:], round(result["summary"]["mean_fst"], 4))
```

- **示例输出：**

```text
2600 ['freq_PopA(n=100)', 'freq_PopB(n=100)', 'Fst'] -0.0001
```

## 27) `dnakit.gwas.gcta_simulate_phenotypes()`

- **作用：** 由真实基因型和指定的因果位点模拟表型。
- **处理方法：** `trait='quantitative'` 运行 `--simu-qt`，`trait='case_control'` 运行 `--simu-cc`。每次运行重新抽随机数，结果不可重复。
- **API：** `dnakit.gwas.gcta_simulate_phenotypes(bfile, out_prefix, *, trait='quantitative', n_cases=None, n_controls=None, causal_loci_file=None, heritability=None, prevalence=None, ..., verbose=True)`。
- **输入：** 基因型、新的输出前缀、性状类型和因果位点文件。
- **输出：** 包含 `status`、`backend`、`version`、`command`、`return_code`、`outputs`、`log_tail` 的字典，另有`causal_variants` 与 `n_individuals`；表型在 `<前缀>.phen`。
- **示例代码：**

```python
import tempfile
from pathlib import Path

out = Path(tempfile.mkdtemp())
from dnakit.gwas.gcta_p08_cli_functions import gcta_simulate_phenotypes

q = gcta_simulate_phenotypes("toy", str(out / "q"), causal_loci_file="causal.snplist", heritability=0.5, replicates=2, verbose=False)
print(q["n_individuals"], q["causal_variants"])
cc = gcta_simulate_phenotypes("toy", str(out / "cc"), trait="case_control", n_cases=40, n_controls=60, prevalence=0.3, causal_loci_file="causal.snplist", verbose=False)
values = [line.split()[2] for line in Path(cc["outputs"]["phen"]).read_text().splitlines()]
print(values.count("2"), values.count("1"))
```

- **示例输出：**

```text
200 [{'QTL': 'rs5', 'RefAllele': 'A', 'Frequency': 0.7925, 'Effect': 0.1}, {'QTL': 'rs100', 'RefAllele': 'A', 'Frequency': 0.62, 'Effect': -0.2}]
40 60
```
