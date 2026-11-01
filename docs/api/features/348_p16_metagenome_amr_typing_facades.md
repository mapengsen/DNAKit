# 宏基因组分类、耐药筛查与微生物分型封装

这一组函数位于 `dnakit.microbiome`，封装 [Kaiju](https://github.com/bioinformatics-centre/kaiju)、[Bracken](https://github.com/jenniferlu717/Bracken)、[geNomad](https://github.com/apcamargo/genomad)、[AMRFinderPlus](https://github.com/ncbi/amr)、[ABRicate](https://github.com/tseemann/abricate)、[mlst](https://github.com/tseemann/mlst)、[fastANI](https://github.com/ParBLiSS/FastANI)、[eggNOG-mapper](https://github.com/eggnogdb/eggnog-mapper)、[DIAMOND](https://github.com/bbuchfink/diamond) 以及 Prokka 的两个数据库转换脚本。每个函数以固定的参数白名单、不经 shell 运行一个程序，并返回可 JSON 序列化的字典。Bracken 表合并与 fastANI 绘图用 Python 实现，不需要外部程序。

DNAKit 不安装这些程序；请把它们放到 `PATH` 或传入 `*_path` 参数，缺少程序时抛出带安装提示的 `BackendUnavailableError`。输出写入新路径，不会覆盖已有文件。下面的示例使用 Kaiju 1.10.3、DIAMOND 2.2.8、fastANI 1.34 与 Bracken 3.1。geNomad、eggNOG-mapper 和 ABRicate 数据库下载需要大型数据库或网络，示例中不运行；它们的命令组装由使用记录型替身程序的单元测试覆盖。

## 示例数据

示例运行前先写入这里声明的输入文件；其余示例在代码中自行生成数据。

**示例文件：** `s1.bracken`

```text
name	taxonomy_id	taxonomy_lvl	kraken_assigned_reads	added_reads	new_est_reads	fraction_total_reads
Species_one	100	S	30	0	30	0.5
Species_two	200	S	30	0	30	0.5
```

**示例文件：** `s2.bracken`

```text
name	taxonomy_id	taxonomy_lvl	kraken_assigned_reads	added_reads	new_est_reads	fraction_total_reads
Species_one	100	S	10	0	10	1.0
```

**示例文件：** `nodes.dmp`

```text
1	|	1	|	no rank	|
10	|	1	|	superkingdom	|
60	|	10	|	genus	|
100	|	60	|	species	|
200	|	60	|	species	|
```

**示例文件：** `names.dmp`

```text
1	|	root	|		|	scientific name	|
10	|	Bacteria	|		|	scientific name	|
60	|	Genus_A	|		|	scientific name	|
100	|	Species_one	|		|	scientific name	|
200	|	Species_two	|		|	scientific name	|
```

## 1) `dnakit.microbiome.p16_microbial_cli_functions.add_kaiju_taxon_names()`

- **作用：** 给 Kaiju 逐 read 输出追加 taxon 名称（`kaiju-addTaxonNames`）。
- **处理方法：** 按 `names.dmp` 追加名称，可丢弃未分类 reads、输出完整谱系或指定层级。
- **API：** `dnakit.microbiome.p16_microbial_cli_functions.add_kaiju_taxon_names(kaiju_output, nodes, names, output_path, *, drop_unclassified=False, full_path=False, path_ranks=None, kaiju_names_path=None, verbose=True)`。
- **输入：** Kaiju 输出、`nodes.dmp`、`names.dmp`、新输出路径。
- **输出：** 字典，含 `status`、`backend`、`version`、`command`、`return_code`、`outputs`（写出的文件）、`log_tail` 和 `provenance`；另有 `line_count`。
- **示例见：** `summarize_kaiju_table()`

## 2) `dnakit.microbiome.p16_microbial_cli_functions.annotate_dna_eggnog_mapper()`

- **作用：** 用 eggNOG-mapper 对 DNA 做直系同源注释（`emapper.py`）。
- **处理方法：** CDS、基因组或宏基因组输入；基因组类先用 prodigal 或同源搜索预测基因，可把注释写回 GFF。
- **API：** `dnakit.microbiome.p16_microbial_cli_functions.annotate_dna_eggnog_mapper(fasta, output_directory, *, prefix='emapper', input_type='CDS', gene_prediction=None, translate=False, decorate_gff=None, data_directory=None, threads=1, emapper_path=None, verbose=True)`。
- **输入：** FASTA、新输出目录、eggNOG 数据目录。
- **输出：** 字典，含 `status`、`backend`、`version`、`command`、`return_code`、`outputs`（写出的文件）、`log_tail` 和 `provenance`；另有 `input_type`、`files`、`rows`。

## 3) `dnakit.microbiome.p16_microbial_cli_functions.build_bracken_kmer_distribution()`

- **作用：** 为 Kraken 数据库构建读长 k-mer 分布（`bracken-build`）。
- **处理方法：** 依次运行 `kmer2read_distr` 模拟读长 k-mer 的分类，再用 `generate_kmer_distribution.py` 汇总；文件写入数据库目录，已存在时拒绝覆盖。
- **API：** `dnakit.microbiome.p16_microbial_cli_functions.build_bracken_kmer_distribution(database, *, read_length=100, kmer_length=35, threads=1, kraken_type='kraken2', bracken_build_path=None, verbose=True)`。
- **输入：** Kraken 数据库目录、读长与 k-mer 长度。
- **输出：** 字典，含 `status`、`backend`、`version`、`command`、`return_code`、`outputs`（写出的文件）、`log_tail` 和 `provenance`；另有 `read_length`；输出 `database<L>mers.kraken` 与 `.kmer_distrib`。

## 4) `dnakit.microbiome.p16_microbial_cli_functions.check_gff_annotation_amrfinder()`

- **作用：** 核对基因注释 GFF 与核酸、蛋白 FASTA 是否一致（AMRFinderPlus `gff_check`）。
- **处理方法：** 按注释格式检查 GFF 的 CDS/gene 与蛋白 ID、位置；退出码 0 表示一致。
- **API：** `dnakit.microbiome.p16_microbial_cli_functions.check_gff_annotation_amrfinder(gff, protein_fasta, nucleotide_fasta, *, annotation_format='genbank', gff_check_path=None, verbose=True)`。
- **输入：** GFF、蛋白 FASTA、核酸 FASTA 与注释格式。
- **输出：** 字典，含 `status`、`backend`、`version`、`command`、`return_code`、`outputs`（写出的文件）、`log_tail` 和 `provenance`；另有 `gff`。

## 5) `dnakit.microbiome.p16_microbial_cli_functions.check_nucleotide_fasta_amrfinder()`

- **作用：** 验证核酸 FASTA 并报告序列数与长度（AMRFinderPlus `fasta_check`）。
- **处理方法：** 程序输出序列数和累计长度，函数解析为 `sequence_count`、`lengths`、`cumulative_lengths`；格式错误时报错。
- **API：** `dnakit.microbiome.p16_microbial_cli_functions.check_nucleotide_fasta_amrfinder(fasta, *, fasta_check_path=None, work_directory=None, verbose=True)`。
- **输入：** 核酸 FASTA。
- **输出：** 字典，含 `status`、`backend`、`version`、`command`、`return_code`、`outputs`（写出的文件）、`log_tail` 和 `provenance`；另有 `sequence_count`、`lengths`、`cumulative_lengths`。

## 6) `dnakit.microbiome.p16_microbial_cli_functions.combine_bracken_tables()`

- **作用：** 把各样本的 Bracken 丰度表合并成多样本表（`combine_bracken_outputs.py`）。
- **处理方法：** 纯 Python 实现，语义与原脚本一致：要求层级相同、taxon ID 一致，输出每个样本的 reads 数和占比（五位小数）。
- **API：** `dnakit.microbiome.p16_microbial_cli_functions.combine_bracken_tables(tables, output_path, *, names=None, verbose=True)`。
- **输入：** Bracken 输出表列表（不是报告）和新输出路径。
- **输出：** 字典，含 `status`、`backend`、`version`、`command`、`return_code`、`outputs`（写出的文件）、`log_tail` 和 `provenance`；另有 `samples`、`taxa`。
- **示例代码：**

```python
from dnakit.microbiome.p16_microbial_cli_functions import combine_bracken_tables

result = combine_bracken_tables(["s1.bracken", "s2.bracken"], "merged.tsv", names=["A", "B"], verbose=False)
print(result["samples"], result["taxa"])
print(open("merged.tsv").read())
```

- **示例输出：**

```text
['A', 'B'] 2
name	taxonomy_id	taxonomy_lvl	A_num	A_frac	B_num	B_frac
Species_one	100	S	30	0.50000	10	1.00000
Species_two	200	S	30	0.50000	0	0.00000
```

## 7) `dnakit.microbiome.p16_microbial_cli_functions.convert_kaiju_to_krona()`

- **作用：** 把 Kaiju 输出转成 Krona 文本（`kaiju2krona`）。
- **处理方法：** 用 `nodes.dmp` 与 `names.dmp` 生成每个 taxon 的路径与计数，可限定层级并计入未分类 reads。
- **API：** `dnakit.microbiome.p16_microbial_cli_functions.convert_kaiju_to_krona(kaiju_output, nodes, names, output_path, *, ranks=None, include_unclassified=False, kaiju2krona_path=None, verbose=True)`。
- **输入：** Kaiju 输出、`nodes.dmp`、`names.dmp`、新输出路径。
- **输出：** 字典，含 `status`、`backend`、`version`、`command`、`return_code`、`outputs`（写出的文件）、`log_tail` 和 `provenance`；另有 `line_count`。
- **示例见：** `summarize_kaiju_table()`

## 8) `dnakit.microbiome.p16_microbial_cli_functions.detect_amr_genes_nucleotide()`

- **作用：** 检测核酸 contig 中的耐药、应激和毒力基因及点突变（`amrfinder --nucleotide`）。
- **处理方法：** 不指定 `organism` 时报告基因；指定后还检测该物种特异的点突变，可加 `--plus` 基因与自定义一致性/覆盖度阈值。
- **API：** `dnakit.microbiome.p16_microbial_cli_functions.detect_amr_genes_nucleotide(contigs, output_path, *, organism=None, database=None, plus=False, min_identity=None, min_coverage=None, translation_table=11, name=None, threads=1, amrfinder_path=None, verbose=True)`。
- **输入：** contig FASTA、新输出路径；AMRFinderPlus 数据库需已安装。
- **输出：** 字典，含 `status`、`backend`、`version`、`command`、`return_code`、`outputs`（写出的文件）、`log_tail` 和 `provenance`；另有 `organism`、`rows`（报告表）。

## 9) `dnakit.microbiome.p16_microbial_cli_functions.download_abricate_database()`

- **作用：** 下载并格式化 ABRicate 参考数据库（`abricate-get_db`，需要网络）。
- **处理方法：** 从公开来源取得所选数据库并用 BLAST+ 建库，写入数据库目录。
- **API：** `dnakit.microbiome.p16_microbial_cli_functions.download_abricate_database(database, *, database_directory=None, force=False, abricate_get_db_path=None, verbose=True)`。
- **输入：** 数据库名，可选目录。
- **输出：** 字典，含 `status`、`backend`、`version`、`command`、`return_code`、`outputs`（写出的文件）、`log_tail` 和 `provenance`；另有 `database`。

## 10) `dnakit.microbiome.p16_microbial_cli_functions.genbank_to_protein_fasta_prokka()`

- **作用：** 把 GenBank 的 CDS 特征转成蛋白 FASTA 参考（Prokka `prokka-genbank_to_fasta_db`）。
- **处理方法：** 翻译每个 CDS 并写成 `>id EC~~~gene~~~product`，默认排除假定蛋白和假基因。
- **API：** `dnakit.microbiome.p16_microbial_cli_functions.genbank_to_protein_fasta_prokka(genbank_files, output_path, *, id_tag=None, include_hypothetical=False, include_pseudo=False, min_length=0, genetic_code=None, separator='~~~', blank='', prokka_script_path=None, verbose=True)`。
- **输入：** GenBank 文件列表与新输出路径。
- **输出：** 字典，含 `status`、`backend`、`version`、`command`、`return_code`、`outputs`（写出的文件）、`log_tail` 和 `provenance`；另有 `proteins`。

## 11) `dnakit.microbiome.p16_microbial_cli_functions.map_fastani_conserved_regions()`

- **作用：** 计算一对基因组的 ANI 并输出双向保守区映射记录（`fastANI --visualize`）。
- **处理方法：** 按片段长度切分查询基因组并映射到参考，摘要给出 ANI 与匹配片段数，`.visual` 文件含每个片段的映射坐标。
- **API：** `dnakit.microbiome.p16_microbial_cli_functions.map_fastani_conserved_regions(query, reference, output_path, *, kmer=16, fragment_length=3000, min_fraction=0.2, threads=1, fastani_path=None, verbose=True)`。
- **输入：** 查询与参考基因组 FASTA、新输出路径。
- **输出：** 字典，含 `status`、`backend`、`version`、`command`、`return_code`、`outputs`（写出的文件）、`log_tail` 和 `provenance`；另有 `ani`、`matched_fragments`、`total_fragments`、`visual_file`。
- **示例代码：**

```python
import random
from dnakit.microbiome.p16_microbial_cli_functions import map_fastani_conserved_regions, plot_fastani_conserved_regions

rng = random.Random(3)
genome = "".join(rng.choice("ACGT") for _ in range(12000))
open("a.fa", "w").write(">a\n" + genome + "\n")
mutated = list(genome)
for position in range(500, 12000, 400):
    mutated[position] = "A" if mutated[position] != "A" else "C"
open("b.fa", "w").write(">b\n" + "".join(mutated) + "\n")
result = map_fastani_conserved_regions("b.fa", "a.fa", "ani.txt", verbose=False)
print(round(result["ani"], 2), result["matched_fragments"], result["total_fragments"])
plot = plot_fastani_conserved_regions(result["visual_file"], "ani.png", verbose=False)
print(plot["fragments"], open("ani.png", "rb").read(4) == b"\x89PNG")
```

- **示例输出：**

```text
99.66 4 4
4 True
```

## 12) `dnakit.microbiome.p16_microbial_cli_functions.merge_kaiju_outputs()`

- **作用：** 合并两个 Kaiju/Kraken 输出的逐 read 分类（`kaiju-mergeOutputs`）。
- **处理方法：** 两个文件按 read 名排序；冲突按 `1`、`2`、`lca` 或 `lowest` 解决，可用得分决定优先。
- **API：** `dnakit.microbiome.p16_microbial_cli_functions.merge_kaiju_outputs(first, second, output_path, *, conflict='lca', nodes=None, use_score=False, kaiju_merge_path=None, verbose=True)`。
- **输入：** 两个排序后的输出、新输出路径，`lca` 时还需 `nodes.dmp`。
- **输出：** 字典，含 `status`、`backend`、`version`、`command`、`return_code`、`outputs`（写出的文件）、`log_tail` 和 `provenance`；另有 `conflict`、`rows`。
- **示例见：** `summarize_kaiju_table()`

## 13) `dnakit.microbiome.p16_microbial_cli_functions.plot_fastani_conserved_regions()`

- **作用：** 绘制两个基因组之间的保守区点图（来自 fastANI `.visual` 文件）。
- **处理方法：** 不运行原 R 脚本，用 matplotlib 画出每个映射片段在参考与查询上的区间。
- **API：** `dnakit.microbiome.p16_microbial_cli_functions.plot_fastani_conserved_regions(visual_file, image_path, *, min_identity=0.0, verbose=True)`。
- **输入：** `.visual` 文件与新的图片路径。
- **输出：** 字典，含 `status`、`backend`、`version`、`command`、`return_code`、`outputs`（写出的文件）、`log_tail` 和 `provenance`；另有 `image`、`fragments`、`rows`。
- **示例见：** `map_fastani_conserved_regions()`

## 14) `dnakit.microbiome.p16_microbial_cli_functions.reestimate_bracken_report()`

- **作用：** 用 Bracken 重新估计丰度并写出带重估 reads 数的 Kraken 风格报告（`bracken -w`）。
- **处理方法：** 对 Kraken 2 报告按层级重新分配 reads，同时写出丰度表和新报告。
- **API：** `dnakit.microbiome.p16_microbial_cli_functions.reestimate_bracken_report(kraken_report, database, output_table, output_report, *, read_length=100, level='S', threshold=0, bracken_path=None, verbose=True)`。
- **输入：** Kraken 2 报告、含 `database<L>mers.kmer_distrib` 的数据库目录、两个新输出路径。
- **输出：** 字典，含 `status`、`backend`、`version`、`command`、`return_code`、`outputs`（写出的文件）、`log_tail` 和 `provenance`；另有 `rows`（name、taxonomy_id、new_est_reads、fraction_total_reads 等）。

## 15) `dnakit.microbiome.p16_microbial_cli_functions.run_genomad()`

- **作用：** 用 geNomad 分类病毒与质粒 contig，可整体运行或分阶段运行。
- **处理方法：** `step` 选择 `end-to-end`、`annotate`、`find-proviruses`、`marker-classification`、`nn-classification`、`aggregated-classification`、`score-calibration` 或 `summary`；各阶段读取同一输出目录中前一阶段的结果。
- **API：** `dnakit.microbiome.p16_microbial_cli_functions.run_genomad(contigs, output_directory, *, step='end-to-end', database=None, threads=1, cleanup=True, splits=None, min_score=None, enable_score_calibration=False, genomad_path=None, verbose=True)`。
- **输入：** contig FASTA、输出目录；需要数据库的阶段还要 geNomad 数据库目录。
- **输出：** 字典，含 `status`、`backend`、`version`、`command`、`return_code`、`outputs`（写出的文件）、`log_tail` 和 `provenance`；另有 `step`、`output_directory`、`files`。

## 16) `dnakit.microbiome.p16_microbial_cli_functions.screen_contigs_abricate()`

- **作用：** 用 BLASTn 在 contig 中筛查耐药、毒力或质粒基因（`abricate`）。
- **处理方法：** 对所选数据库做 BLASTn，按最小一致性与覆盖度过滤，报告写到文件并解析前 1,000 行。
- **API：** `dnakit.microbiome.p16_microbial_cli_functions.screen_contigs_abricate(contigs, output_path, *, database='resfinder', min_identity=75.0, min_coverage=0.0, data_directory=None, abricate_path=None, verbose=True)`。
- **输入：** contig FASTA、新输出路径、已安装的数据库名。
- **输出：** 字典，含 `status`、`backend`、`version`、`command`、`return_code`、`outputs`（写出的文件）、`log_tail` 和 `provenance`；另有 `database`、`rows`。

## 17) `dnakit.microbiome.p16_microbial_cli_functions.search_dna_translated_diamond()`

- **作用：** 用 DIAMOND blastx 把 DNA 查询翻译后比对蛋白数据库。
- **处理方法：** 支持表格输出、每条 read 的 LCA 分类（输出格式 102）、移码感知、长读长和穷尽搜索 `--swipe`。
- **API：** `dnakit.microbiome.p16_microbial_cli_functions.search_dna_translated_diamond(queries, database, output_path, *, output_format='6', taxonomic_classification=False, frameshift_penalty=None, long_reads=False, exhaustive=False, max_target_seqs=None, evalue=0.001, sensitivity=None, threads=1, diamond_path=None, verbose=True)`。
- **输入：** DNA 查询、`.dmnd` 数据库、新输出路径。
- **输出：** 字典，含 `status`、`backend`、`version`、`command`、`return_code`、`outputs`（写出的文件）、`log_tail` 和 `provenance`；另有 `output_format`、`rows`。
- **示例代码：**

```python
import random, subprocess
codon = {"A": "GCT", "C": "TGT", "D": "GAT", "E": "GAA", "F": "TTT", "G": "GGT", "H": "CAT", "I": "ATT", "K": "AAA", "L": "CTG",
         "M": "ATG", "N": "AAT", "P": "CCT", "Q": "CAA", "R": "CGT", "S": "TCT", "T": "ACT", "V": "GTT", "W": "TGG", "Y": "TAT"}
rng = random.Random(5)
proteins = {f"p{i}_{taxon}": "".join(rng.choice(list(codon)) for _ in range(300)) for i, taxon in enumerate((100, 200))}
open("db.faa", "w").write("".join(f">{k}\n{v}\n" for k, v in proteins.items()))
reads = {f"read_{k}_{j}": "".join(codon[a] for a in v[j * 30:j * 30 + 50]) for k, v in proteins.items() for j in range(3)}
reads["read_random"] = "".join(rng.choice("ACGT") for _ in range(150))
open("reads.fa", "w").write("".join(f">{k}\n{v}\n" for k, v in reads.items()))
import subprocess
from dnakit.microbiome.p16_microbial_cli_functions import search_dna_translated_diamond

subprocess.run(["diamond", "makedb", "--in", "db.faa", "-d", "proteins", "--quiet"], check=True, capture_output=True)
result = search_dna_translated_diamond("reads.fa", "proteins.dmnd", "hits.tsv", exhaustive=True, verbose=False)
for row in result["rows"][:3]:
    print(row["column_1"], row["column_2"], row["column_3"])
```

- **示例输出：**

```text
read_p0_100_0 p0_100 100
read_p0_100_1 p0_100 100
read_p0_100_2 p0_100 100
```

## 18) `dnakit.microbiome.p16_microbial_cli_functions.search_kaiju_best_match()`

- **作用：** 为 DNA reads 查找最匹配的数据库蛋白，不做分类（`kaijux`）。
- **处理方法：** 用 Kaiju 的 `.fmi` 蛋白索引做翻译搜索，输出 C/U、read 名、得分和匹配的数据库标识。
- **API：** `dnakit.microbiome.p16_microbial_cli_functions.search_kaiju_best_match(reads, database, output_path, *, reads_2=None, mode='greedy', mismatches=3, min_match_length=11, min_score=65, max_evalue=0.01, seg_filter=True, threads=1, kaijux_path=None, verbose=True)`。
- **输入：** reads（FASTA/FASTQ）、`.fmi` 索引、新输出路径。
- **输出：** 字典，含 `status`、`backend`、`version`、`command`、`return_code`、`outputs`（写出的文件）、`log_tail` 和 `provenance`；另有 `rows`（classified、read、score、matches）、`classified`。
- **示例代码：**

```python
import random, subprocess
codon = {"A": "GCT", "C": "TGT", "D": "GAT", "E": "GAA", "F": "TTT", "G": "GGT", "H": "CAT", "I": "ATT", "K": "AAA", "L": "CTG",
         "M": "ATG", "N": "AAT", "P": "CCT", "Q": "CAA", "R": "CGT", "S": "TCT", "T": "ACT", "V": "GTT", "W": "TGG", "Y": "TAT"}
rng = random.Random(5)
proteins = {f"p{i}_{taxon}": "".join(rng.choice(list(codon)) for _ in range(300)) for i, taxon in enumerate((100, 200))}
open("db.faa", "w").write("".join(f">{k}\n{v}\n" for k, v in proteins.items()))
reads = {f"read_{k}_{j}": "".join(codon[a] for a in v[j * 30:j * 30 + 50]) for k, v in proteins.items() for j in range(3)}
reads["read_random"] = "".join(rng.choice("ACGT") for _ in range(150))
open("reads.fa", "w").write("".join(f">{k}\n{v}\n" for k, v in reads.items()))
subprocess.run(["kaiju-mkbwt", "-n", "2", "-l", "1", "-a", "ACDEFGHIKLMNPQRSTVWY", "-o", "db", "db.faa"], check=True, capture_output=True)
subprocess.run(["kaiju-mkfmi", "db"], check=True, capture_output=True)
from dnakit.microbiome.p16_microbial_cli_functions import search_kaiju_best_match

result = search_kaiju_best_match("reads.fa", "db.fmi", "best.out", verbose=False)
print(result["classified"])
for row in result["rows"][:3] + result["rows"][-1:]:
    print(f'{row["classified"]} {row["read"]} {row.get("matches", "")}'.rstrip())
```

- **示例输出：**

```text
6
C read_p0_100_0 p0_100,
C read_p0_100_1 p0_100,
C read_p0_100_2 p0_100,
U read_random
```

## 19) `dnakit.microbiome.p16_microbial_cli_functions.summarize_abricate_reports()`

- **作用：** 把多个 ABRicate 报告汇总成基因有/无矩阵（`abricate --summary`）。
- **处理方法：** 合并报告，列出每个样本中每个基因的覆盖度。
- **API：** `dnakit.microbiome.p16_microbial_cli_functions.summarize_abricate_reports(reports, output_path, *, abricate_path=None, verbose=True)`。
- **输入：** ABRicate 报告列表、新输出路径。
- **输出：** 字典，含 `status`、`backend`、`version`、`command`、`return_code`、`outputs`（写出的文件）、`log_tail` 和 `provenance`；另有 `rows`。

## 20) `dnakit.microbiome.p16_microbial_cli_functions.summarize_kaiju_table()`

- **作用：** 按分类层级汇总各样本的 Kaiju 分类（`kaiju2table`）。
- **处理方法：** 统计每个 taxon 的 reads 与百分比，可设最小百分比或 reads 数、展开病毒、输出完整路径。
- **API：** `dnakit.microbiome.p16_microbial_cli_functions.summarize_kaiju_table(kaiju_outputs, nodes, names, output_path, *, rank='species', min_percent=None, min_reads=None, expand_viruses=False, exclude_unclassified=False, full_path=False, path_ranks=None, kaiju2table_path=None, verbose=True)`。
- **输入：** 一个或多个 Kaiju 输出、`nodes.dmp`、`names.dmp`、层级。
- **输出：** 字典，含 `status`、`backend`、`version`、`command`、`return_code`、`outputs`（写出的文件）、`log_tail` 和 `provenance`；另有 `rank`、`rows`（file、percent、reads、taxon_id、taxon_name）。
- **示例代码：**

```python
import subprocess
from dnakit.microbiome.p16_microbial_cli_functions import add_kaiju_taxon_names, merge_kaiju_outputs, summarize_kaiju_table

subprocess.run(["kaiju", "-t", "nodes.dmp", "-f", "db.fmi", "-i", "reads.fa", "-o", "kaiju.out"], check=True, capture_output=True)
table = summarize_kaiju_table(["kaiju.out"], "nodes.dmp", "names.dmp", "species.tsv", rank="species", verbose=False)
for row in table["rows"]:
    if row["reads"]:
        print(row["taxon_name"], row["reads"])
add_kaiju_taxon_names("kaiju.out", "nodes.dmp", "names.dmp", "named.out", verbose=False)
print(open("named.out").read().splitlines()[0].split("\t")[-1])
merged = merge_kaiju_outputs("kaiju.out", "kaiju.out", "merged.out", nodes="nodes.dmp", verbose=False)
print(open("merged.out").read() == open("kaiju.out").read())
```

- **示例输出：**

```text
Species_one 3
Species_two 3
unclassified 1
Species_one
True
```

## 21) `dnakit.microbiome.p16_microbial_cli_functions.translate_abricate_genes_prokka()`

- **作用：** 把 ABRicate 核酸基因 FASTA 翻译成 Prokka 蛋白参考（`prokka-abricate_to_fasta_db`）。
- **处理方法：** 强制 Met 起始、去掉终止密码子；序列含内部终止或蛋白重复时程序报错。
- **API：** `dnakit.microbiome.p16_microbial_cli_functions.translate_abricate_genes_prokka(abricate_sequences, output_path, *, prokka_script_path=None, verbose=True)`。
- **输入：** ABRicate 序列文件与新输出路径。
- **输出：** 字典，含 `status`、`backend`、`version`、`command`、`return_code`、`outputs`（写出的文件）、`log_tail` 和 `provenance`；另有 `proteins`。

## 22) `dnakit.microbiome.p16_microbial_cli_functions.type_sequence_type_mlst()`

- **作用：** 对组装 contig 做多位点序列分型（`mlst`）。
- **处理方法：** 自动检测方案或强制指定，输出 ST 与各位点等位基因；新等位基因可导出为 FASTA。
- **API：** `dnakit.microbiome.p16_microbial_cli_functions.type_sequence_type_mlst(contigs, output_path, *, scheme=None, novel_alleles_path=None, min_identity=95.0, min_coverage=10.0, min_score=50.0, threads=1, mlst_path=None, verbose=True)`。
- **输入：** contig FASTA、新输出路径。
- **输出：** 字典，含 `status`、`backend`、`version`、`command`、`return_code`、`outputs`（写出的文件）、`log_tail` 和 `provenance`；另有 `fields`（file、scheme、ST、等位基因）。

## 23) `dnakit.microbiome.p16_microbial_cli_functions.type_shiga_toxin_operons()`

- **作用：** 确定组装序列中的志贺毒素 stx 操纵子类型（`stxtyper`）。
- **处理方法：** 对核酸序列比对 stx 参考，输出操纵子完整性与亚型。
- **API：** `dnakit.microbiome.p16_microbial_cli_functions.type_shiga_toxin_operons(contigs, output_path, *, amrfinder_format=False, name=None, threads=1, stxtyper_path=None, verbose=True)`。
- **输入：** 核酸 FASTA、新输出路径。
- **输出：** 字典，含 `status`、`backend`、`version`、`command`、`return_code`、`outputs`（写出的文件）、`log_tail` 和 `provenance`；另有 `rows`（stx_type、operon、一致性）。
