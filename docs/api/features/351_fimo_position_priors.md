# FIMO 位置特异先验的基序扫描

这一组只有一个函数，位于 `dnakit.regulation`：对 [MEME suite](https://meme-suite.org) 中 `fimo`（测试版本 5.5.9）的带类型命令行封装。它在 DNA 序列中扫描 MEME 格式基序的出现位置；给出位置特异先验（PSP）时，用 Cuellar-Partida 等人 2012（Bioinformatics 28:56-62）的对数后验几率代替普通对数几率给每个出现位置打分。先验可以来自 DNase 或组蛋白标记等表观基因组信号，所以开放区域内的基序匹配会比关闭区域内的同一匹配排名更高。

DNAKit 不内置 MEME suite。请自行安装（`conda install -c bioconda meme`，或从 meme-suite.org 下载发布版），把 `fimo` 放到 `PATH`，或传入 `fimo_path`。先验分布文件由同一套件的 `compute-prior-dist` 或 `create-priors` 生成。找不到 `fimo` 时函数会报错并给出安装提示。函数不使用 shell，也不访问网络。

该函数对应工作清单中的“带位置特异先验的基序出现位置扫描”（`--psp` 与 `--prior-dist`）。不带先验的普通 PWM 扫描已由 `dnakit.predictions` 中的 `scan_motifs_with_pwms()` 提供。

## 示例数据

一条 60 bp 的序列，`GATCAG` 位点位于第 13、22、34 位；一个识别 `GATCAG` 的基序；一个只在第 13 位提高先验的 wiggle 先验文件（没有给值的位置使用先验中位数）；以及分箱后的先验分布（最小值、最大值、中位数，然后是四个箱中各自的先验占比）。

**示例文件：** `example_seq.fa`

```text
>seq1
ACGTTGCAACCTGATCAGTAGGATCAGCTTAACGATCAGGTCCATTGACAGTCCGATGCA
```

**示例文件：** `example_motif.meme`

```text
MEME version 4

ALPHABET= ACGT

strands: + -

Background letter frequencies
A 0.25 C 0.25 G 0.25 T 0.25

MOTIF TESTMOTIF
letter-probability matrix: alength= 4 w= 6 nsites= 20 E= 0
0.05 0.05 0.85 0.05
0.85 0.05 0.05 0.05
0.05 0.05 0.05 0.85
0.05 0.85 0.05 0.05
0.85 0.05 0.05 0.05
0.05 0.05 0.85 0.05
```

**示例文件：** `example.wig`

```text
variableStep chrom=seq1 span=1
12 0.05
13 0.05
```

**示例文件：** `example.dist`

```text
#min 0.0100
#max 0.0500
#median 0.0100
0.8
0.0
0.0
0.2
```

## 1) `dnakit.regulation.scan_motifs_with_position_priors()`

- **作用：** 用 FIMO 在 DNA 中查找基序出现位置，并用位置特异先验（对数后验几率）打分，同时给出该分数的 p 值和 q 值。
- **处理方法：** 运行 `fimo --psp <prior_file> --prior-dist <prior_distribution_file> --alpha <prior_alpha> --thresh <p_value_threshold> ... <motif_file> <sequence_file>`，输出到新建文件夹，再解析 `fimo.tsv`。与 FIMO 一致，两个先验文件必须同时提供。默认值与 FIMO 相同：p 值阈值 1e-4、扫描双链、基序伪计数 0.1、先验 alpha 为 1.0、保存 100,000 个分数。可选参数：`use_q_value_threshold`（`--qv-thresh`）、`compute_q_values=False`（`--no-qvalue`）、`scan_reverse_complement=False`（`--norc`）、`max_strand`、`motif_ids`（`--motif`，可重复）、`motif_pseudocount`、`background_file` 或 `background_source`（`motif`、`uniform`、`nrdb`，对应 `--bfile`）、`parse_genomic_coordinates=False`（`--no-pgc`）和 `max_stored_scores`。wiggle 先验的坐标是位点的 1 起始起点。
- **API：** `dnakit.regulation.scan_motifs_with_position_priors(motif_file, sequence_file, prior_file, prior_distribution_file, output_dir, *, prior_alpha=1.0, p_value_threshold=1e-4, use_q_value_threshold=False, compute_q_values=True, scan_reverse_complement=True, max_strand=False, motif_ids=None, motif_pseudocount=0.1, background_file=None, background_source=None, parse_genomic_coordinates=True, max_stored_scores=100000, max_rows=1000, fimo_path=None)`。
- **输入：** MEME 格式的基序文件、FASTA 文件、wiggle 或 MEME PSP 格式的先验、其分箱分布文件，以及一个尚不存在的文件夹。
- **输出：** 字典，含 `status`（`completed` 或 `failed`）、`backend`、`version`、`command`、`return_code`、`log_tail`、`output_dir`、`outputs`（`fimo.tsv`、`fimo.gff`、`fimo.html`、`fimo.xml`、`cisml.xml`、`best_site.narrowPeak`）、`match_count`、`matches`（前 `max_rows` 行，含 `motif_id`、`motif_alt_id`、`sequence_name`、`start`、`stop`、`strand`、`score`、`p_value`、`q_value`、`matched_sequence`）和 `matches_truncated`。
- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.regulation.fimo_p17_cli_functions import scan_motifs_with_position_priors

output = Path(tempfile.mkdtemp()) / "fimo_out"
result = scan_motifs_with_position_priors(
    "example_motif.meme",
    "example_seq.fa",
    "example.wig",
    "example.dist",
    str(output),
    p_value_threshold=1e-3,
    verbose=False,
)
print(result["status"], result["backend"], result["match_count"])
for row in result["matches"]:
    print(row["motif_id"], row["start"], row["stop"], row["strand"], row["score"], row["p_value"], row["matched_sequence"])
print(sorted(Path(path).name for path in result["outputs"]))
```

- **示例输出：**

```text
completed fimo 4
TESTMOTIF 13 18 + 6.35714 4.74e-05 GATCAG
TESTMOTIF 11 16 - 3.98214 0.000237 GATCAG
TESTMOTIF 22 27 + 3.98214 0.000237 GATCAG
TESTMOTIF 34 39 + 3.98214 0.000237 GATCAG
['best_site.narrowPeak', 'cisml.xml', 'fimo.gff', 'fimo.html', 'fimo.tsv', 'fimo.xml']
```
