# 系统发育与分子进化：参数与补充说明

用法和示例见[功能教程](../../api/features/30_phylogeny.md)。

各功能按作用、API、输入、输出、示例代码和结果说明介绍。运行前按说明准备凭据、本地文件与可选依赖；远程结果随数据库版本与查询时间变化。

## 计算两条编码序列的 dN/dS {#dnakit.databases.sequence_dn_ds-parameters}

`dnakit.databases.sequence_dn_ds(*, seq1=None, seq2=None, fasta1_path=None, fasta2_path=None, verbose=True)`

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `seq1` | string &#124; null | None | First coding sequence (in-frame, codon-aligned to seq2). DNA/RNA letters. |
| `seq2` | string &#124; null | None | Second coding sequence (same length/frame as seq1). |
| `fasta1_path` | string &#124; null | None | Alternative to seq1: path to a single-record FASTA |
| `fasta2_path` | string &#124; null | None | Alternative to seq2: path to a single-record FASTA |
| `verbose` | bool | True | 运行提示开关。 |
