# CRISPR 筛选分析、脱靶搜索与 DNABERT 嵌入

`dnakit.design_crispr` 新增 [MAGeCK](https://sourceforge.net/p/mageck/wiki/Home/)（count、RRA 检验、MLE、通路、绘图）、[DrugZ](https://github.com/hart-lab/drugz) 和 [Cas-OFFinder](https://github.com/snugel/cas-offinder) 的封装；`dnakit.representations` 新增 `embed_dna_with_dnabert`，复现 LucaOneTasks 的 DNABERT-S / DNABERT-2 推理脚本（token 矩阵或 `[CLS]` 向量）。DNAKit 不安装这些程序；请放到 `PATH` 或传入 `*_path`。DrugZ 是脚本：下载 `drugz.py` 并作为 `drugz_script` 传入。测试版本：MAGeCK 0.5.9.5、Cas-OFFinder 2.4.1（通过 pocl 使用 CPU）、DrugZ 1.1.0.2。DNABERT 函数会加载检查点自带的代码，需要 `allow_remote_code=True` 和本地检查点；它与 LucaOneTasks 脚本对照，矩阵完全一致（最大差 0）。

## 1) `dnakit.design_crispr.p16_crispr_screen_cli_functions.count_sgrnas_mageck()`

- **作用：** 对 sgRNA 文库统计 FASTQ（或比对）中的 sgRNA 读数（`mageck count`）。
- **处理方法：** 输出计数表、归一化计数表与计数摘要。
- **API：** `dnakit.design_crispr.p16_crispr_screen_cli_functions.count_sgrnas_mageck(library, output_prefix, *, fastq=None, sample_labels=None, alignment_files=None, normalization='median', trim_5=None, sgrna_length=None, reverse_complement=False, count_n=False, mageck_path=None, verbose=True)`。
- **输入：** sgRNA 文库、FASTQ 列表、样本标签、输出前缀。
- **输出：** 字典，含 `status`、`backend`、`version`、`command`、`return_code`、`outputs`（写出的文件）、`log_tail` 和 `provenance`；另有 `summary`、`count_table`。
- **示例代码：**

```python
import random
from dnakit.design_crispr.p16_crispr_screen_cli_functions import count_sgrnas_mageck, rank_genes_mageck

rng = random.Random(2)
genes = [f"G{i}" for i in range(40)]
library = [(f"{g}_sg{j}", "".join(rng.choice("ACGT") for _ in range(20)), g) for g in genes for j in range(4)]
open("library.txt", "w").write("".join(f"{a}\t{b}\t{c}\n" for a, b, c in library))
for name, factor in (("c1", 1.0), ("c2", 1.0), ("t1", 0.1), ("t2", 0.1)):
    with open(name + ".fq", "w") as handle:
        n = 0
        for _sid, seq, gene in library:
            count = rng.randint(80, 150)
            for _ in range(int(count * factor) if gene in ("G1", "G2") else count):
                handle.write(f"@r{n}\nTT{seq}AA\n+\n{'I' * 24}\n")
                n += 1
labels = ["c1", "c2", "t1", "t2"]
counted = count_sgrnas_mageck("library.txt", "screen", fastq=[n + ".fq" for n in labels], sample_labels=labels, trim_5=2, verbose=False)
print([(row["Label"], row["Mapped"] == row["Reads"]) for row in counted["summary"]])
ranked = rank_genes_mageck(counted["count_table"], "rra", treatment=["t1", "t2"], control=["c1", "c2"], verbose=False)
for row in ranked["genes"][:2]:
    print(row["id"], row["neg|rank"], row["neg|fdr"])
```

- **示例输出：**

```text
[('c1', True), ('c2', True), ('t1', True), ('t2', True)]
G2 1 0.0001
G1 2 0.0001
```

## 2) `dnakit.design_crispr.p16_crispr_screen_cli_functions.enrich_pathways_mageck()`

- **作用：** 对排序基因做通路/基因集富集（`mageck pathway`，GSEA 或 RRA）。
- **处理方法：** 输入 `test`/`mle` 的基因摘要和 GMT 文件。
- **API：** `dnakit.design_crispr.p16_crispr_screen_cli_functions.enrich_pathways_mageck(gene_ranking, gene_sets, output_prefix, *, method='gsea', single_ranking=False, ranking_column=None, sort_criteria=None, permutations=None, mageck_path=None, verbose=True)`。
- **输入：** 基因排序文件、GMT、输出前缀。
- **输出：** 字典，含 `status`、`backend`、`version`、`command`、`return_code`、`outputs`（写出的文件）、`log_tail` 和 `provenance`；另有 `pathways`。

## 3) `dnakit.design_crispr.p16_crispr_screen_cli_functions.estimate_gene_essentiality_mageck()`

- **作用：** 用设计矩阵估计基因必需性（`mageck mle`，EM 与置换 p 值）。
- **处理方法：** 同样支持拷贝数校正与估计。
- **API：** `dnakit.design_crispr.p16_crispr_screen_cli_functions.estimate_gene_essentiality_mageck(count_table, design_matrix, output_prefix, *, beta_labels=None, include_samples=None, normalization='median', cnv_normalization=None, cell_line=None, cnv_estimation=None, control_sgrna=None, permutation_rounds=None, remove_outliers=False, threads=1, mageck_path=None, verbose=True)`。
- **输入：** 计数表、设计矩阵、输出前缀。
- **输出：** 字典，含 `status`、`backend`、`version`、`command`、`return_code`、`outputs`（写出的文件）、`log_tail` 和 `provenance`；另有 `genes`、`gene_summary`。

## 4) `dnakit.design_crispr.p16_crispr_screen_cli_functions.plot_screen_results_mageck()`

- **作用：** 为所列基因绘制 sgRNA 计数图（`mageck plot`）。
- **处理方法：** MAGeCK 调用 R 生成 PDF；未生成 PDF 时报错。
- **API：** `dnakit.design_crispr.p16_crispr_screen_cli_functions.plot_screen_results_mageck(count_table, gene_summary, output_prefix, *, genes, samples=None, normalization='median', mageck_path=None, verbose=True)`。
- **输入：** 计数表、基因摘要、基因列表、输出前缀。
- **输出：** 字典，含 `status`、`backend`、`version`、`command`、`return_code`、`outputs`（写出的文件）、`log_tail` 和 `provenance`；另有 `genes`。

## 5) `dnakit.design_crispr.p16_crispr_screen_cli_functions.score_drug_gene_interactions_drugz()`

- **作用：** 用 DrugZ 评估配对 sgRNA 计数筛选中的药物-基因互作。
- **处理方法：** 默认按顺序配对比较对照与药物样本，`unpaired` 比较均值倍数变化；输出 normZ、p 值与 FDR。
- **API：** `dnakit.design_crispr.p16_crispr_screen_cli_functions.score_drug_gene_interactions_drugz(count_table, output_path, *, control_samples, drug_samples, drugz_script, unpaired=False, fold_change_path=None, remove_genes=None, pseudocount=5, min_observations=1, half_window_size=500, drugz_python=None, verbose=True)`。
- **输入：** 计数表、对照与药物样本、下载的 `drugz.py` 脚本。
- **输出：** 字典，含 `status`、`backend`、`version`、`command`、`return_code`、`outputs`（写出的文件）、`log_tail` 和 `provenance`；另有 `unpaired`、`genes`。

## 6) `dnakit.design_crispr.p16_crispr_screen_cli_functions.search_offtargets_cas_offinder()`

- **作用：** 在基因组中搜索容忍错配和 PAM 的 CRISPR 脱靶位点（`cas-offinder`）。
- **处理方法：** 写出 Cas-OFFinder 输入文件（基因组、模式、向导序列与最大错配数）并运行，解析位点表。
- **API：** `dnakit.design_crispr.p16_crispr_screen_cli_functions.search_offtargets_cas_offinder(genome, guides, output_path, *, pattern='NNNNNNNNNNNNNNNNNNNNNRG', device='C', cas_offinder_path=None, verbose=True)`。
- **输入：** 基因组 FASTA 或目录、向导序列到最大错配数的映射。
- **输出：** 字典，含 `status`、`backend`、`version`、`command`、`return_code`、`outputs`（写出的文件）、`log_tail` 和 `provenance`；另有 `site_count`、`sites`、`input_file`。
- **示例代码：**

```python
import random
from dnakit.design_crispr.p16_crispr_screen_cli_functions import search_offtargets_cas_offinder

rng = random.Random(6)
chrom = "".join(rng.choice("ACGT") for _ in range(6000))
guide = chrom[1000:1020]
chrom = chrom[:1020] + "TGG" + chrom[1023:]
open("target.fa", "w").write(">chrT\n" + chrom + "\n")
result = search_offtargets_cas_offinder("target.fa", {guide + "NNN": 2}, "hits.txt", pattern="N" * 20 + "NGG", verbose=False)
exact = [site for site in result["sites"] if site["mismatches"] == 0]
print(exact[0]["chrom"], exact[0]["position"], exact[0]["strand"], exact[0]["site"][-3:])
```

- **示例输出：**

```text
chrT 1000 + TGG
```

## 7) `dnakit.design_crispr.p16_crispr_screen_cli_functions.rank_genes_mageck()`

- **作用：** 对计数表做稳健秩聚合基因排序（`mageck test`），含负选与正选。
- **处理方法：** 可配对比较，可提供拷贝数矩阵做校正（`--cnv-norm`）或估计拷贝数（`--cnv-est`）。
- **API：** `dnakit.design_crispr.p16_crispr_screen_cli_functions.rank_genes_mageck(count_table, output_prefix, *, treatment, control=None, day0_label=None, paired=False, normalization='median', adjust_method='fdr', gene_lfc_method=None, cnv_normalization=None, cell_line=None, cnv_estimation=None, control_sgrna=None, remove_zero=None, write_normalized_counts=False, mageck_path=None, verbose=True)`。
- **输入：** 计数表、处理与对照样本标签、输出前缀。
- **输出：** 字典，含 `status`、`backend`、`version`、`command`、`return_code`、`outputs`（写出的文件）、`log_tail` 和 `provenance`；另有 `genes`、`gene_summary`、`sgrna_summary`。
- **示例见：** `count_sgrnas_mageck()`

## 8) `dnakit.representations.lucaonetasks_p16_functions.embed_dna_with_dnabert()`

- **作用：** 用 DNABERT-S 或 DNABERT-2 为 DNA 序列生成嵌入（LucaOneTasks 协议）。
- **处理方法：** 序列转大写并只保留 A-Z，U 读作 T，按 `truncation_length` 截断后由检查点的分词器编码，返回最后一层隐藏状态的矩阵、`[CLS]` 向量或全部 token 的均值。
- **API：** `dnakit.representations.lucaonetasks_p16_functions.embed_dna_with_dnabert(sequences, *, model='dnabert-s', checkpoint_path=None, embedding_type='matrix', include_special_tokens=False, truncation_length=4094, truncate_side='right', output_path=None, device='auto', allow_remote_code=False, verbose=True)`。
- **输入：** 序列列表或 id 到序列的映射、本地检查点目录，并显式 `allow_remote_code=True`。
- **输出：** 字典，含 `status`、`model`、`embedding_dimension`、`sequences`、`embeddings`、`outputs`、`provenance`。
