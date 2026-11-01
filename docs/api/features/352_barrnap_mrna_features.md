# Barrnap 信使 RNA 特征预测

这一组只有一个函数，位于 `dnakit.annotation`：对 [Barrnap](https://github.com/tseemann/barrnap) 1.x（按 1.10.5 编写）`barrnap --mrna` 的带类型命令行封装。除 rRNA 外，Barrnap 1.x 还用 Pyrodigal 预测蛋白编码基因，并以 GFF3 报告 `mRNA`（向上游延伸覆盖核糖体结合位点）、`RBS`（含其基序）、`CDS`（用 DIAMOND 对所选界的 Barrnap 蛋白库比对命名）、`sig_pep`（信号肽，ZiggyPep）和不依赖 Rho 的 `terminator`（TransTermHP）。

DNAKit 不内置 Barrnap。请连同依赖一起安装（`conda install -c bioconda -c conda-forge barrnap`），把 `barrnap` 放到 `PATH`，或传入 `barrnap_path`。找不到时函数报错并给出安装提示。函数不使用 shell，也不访问网络。

输入合计短于约 19 kb 时，Barrnap 会警告并跳过 mRNA 搜索，所以应提供基因组或重叠群，而不是单个基因。Barrnap 的 rRNA、tRNA、ncRNA 搜索是另外的操作；`include_rrna=True` 只是把 Barrnap 默认的 rRNA 搜索加入同一个 GFF3。

## 示例数据

任何不短于 19 kb 的细菌基因组或重叠群 FASTA，例如 `genome.fa`。下面的示例没有声明该文件，因此文档检查脚本不会运行本页；输出只展示真实结果的形状，不是实测值。

## 1) `dnakit.annotation.predict_mrna_features_barrnap()`

- **作用：** 用 Barrnap 在 DNA 中预测蛋白编码基因及其 RBS、信号肽和终止子特征。
- **处理方法：** 运行 `barrnap --quiet --mrna --no-rrna --kingdom <kingdom> --gcode <genetic_code> --evalue <evalue> --threads <threads> <序列文件>`，把标准输出的 GFF3 流式写入新的 `output_gff`，标准错误单独保存，并解析特征。`kingdom` 取 `bac`、`arc` 或 `fun`；`genetic_code='auto'` 时 Barrnap 分别使用密码子表 11、25、1。`add_ids` 对应 `--addids`，`include_rrna` 把 `--no-rrna` 换成 `--rrna`，`database_dir` 对应 `--dbdir`。
- **API：** `dnakit.annotation.predict_mrna_features_barrnap(sequence_files, output_gff, *, kingdom="bac", genetic_code="auto", evalue=0.001, threads=1, add_ids=False, include_rrna=False, database_dir=None, max_features=1000, barrnap_path=None)`。
- **输入：** 一个或多个 DNA 文件（FASTA、GenBank 或 FASTQ，可 gzip 压缩），以及一个尚不存在的 GFF3 路径。
- **输出：** 字典，含 `status`、`backend`、`version`、`command`、`return_code`、`log_tail`、`gff_file`、`feature_counts`（类型到数量）、`feature_count`、`features`（前 `max_features` 行，含 `seqid`、`source`、`type`、`start`、`end`、`score`、`strand`、`phase` 和 `attributes` 字典）与 `features_truncated`。
- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.annotation.barrnap_p17_cli_functions import predict_mrna_features_barrnap

output = Path(tempfile.mkdtemp()) / "mrna.gff3"
result = predict_mrna_features_barrnap(["genome.fa"], str(output), kingdom="bac", threads=4)
print(result["status"], result["feature_counts"])
```

- **示例输出（仅示意形状，非实测）：**

```text
completed {'mRNA': 190, 'RBS': 170, 'CDS': 190, 'sig_pep': 12, 'terminator': 40}
```
