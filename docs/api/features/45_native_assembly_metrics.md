# 原生 Assembly 连续性指标

这些函数只对用户提供的 contig 长度计算连续性指标，不运行组装器，也不读取参考基因组或 BUSCO 数据库。方法定义可与 Assemblathon 2 的组装评估论文对照：<https://doi.org/10.1038/nmeth.1935>。

## 1) `dnakit.assembly.calculate_n50()`

- **作用：** 计算 N50 和对应的 L50。
- **处理方法：** 将 contig 按长度降序排列，累计长度达到总长度 50% 时取当前 contig。
- **API：** `dnakit.assembly.calculate_n50(input, options=None)`。
- **输入：** `[100, 60, 40]` 或 `{"lengths": [100, 60, 40]}`。
- **输出：** 带 `status`、`result` 和 `provenance` 的 JSON 可序列化字典。
- **示例代码：** `dnakit.assembly.calculate_n50([100, 60, 40])`。
- **示例结果：** `{"result": {"n50": 100, "l50": 1, "total_length": 200}}`。

## 2) `dnakit.assembly.calculate_nx()`

- **作用：** 计算任意百分位的 Nx。
- **处理方法：** 用 `options.x` 或 `options.percentile` 指定累计阈值。
- **API：** `dnakit.assembly.calculate_nx(input, options={"x": 90})`。
- **输入：** 正整数 contig 长度列表。
- **输出：** 包含 `nx`、阈值排名和总长度的结果字典。
- **示例代码：** `dnakit.assembly.calculate_nx([100, 60, 40], options={"x": 90})`。
- **示例结果：** `{"result": {"percentile": 90.0, "nx": 40, "contigs_to_threshold": 3}}`。

## 3) `dnakit.assembly.calculate_l50()`

- **作用：** 计算达到指定累计百分位所需的 contig 数量。
- **处理方法：** 与 Nx 使用同一排序和累计规则，默认百分位为 50。
- **API：** `dnakit.assembly.calculate_l50(input, options=None)`。
- **输入：** 正整数 contig 长度列表。
- **输出：** `l50` 或指定百分位的累计排名。
- **示例代码：** `dnakit.assembly.calculate_l50([100, 60, 40])`。
- **示例结果：** `{"result": {"l50": 1, "n_at_percentile": 100}}`。

## 4) `dnakit.assembly.calculate_ng50()`

- **作用：** 使用预期基因组大小计算 NG50。
- **处理方法：** 用 `options.genome_size` 替代实际 contig 总长度作为累计分母。
- **API：** `dnakit.assembly.calculate_ng50(input, options={"genome_size": 250})`。
- **输入：** contig 长度列表和明确的预期基因组大小。
- **输出：** 包含 `ng50`、预期基因组大小和组装总长度的结果字典。
- **示例代码：** `dnakit.assembly.calculate_ng50([100, 60, 40], options={"genome_size": 250})`。
- **示例结果：** `{"result": {"ng50": 60, "expected_genome_size": 250}}`。

## 5) `dnakit.assembly.calculate_stats()`

- **作用：** 一次返回常用组装长度、N50/L50 和 N90/L90 指标。
- **处理方法：** 在同一组显式 contig 长度上计算描述性统计。
- **API：** `dnakit.assembly.calculate_stats(input, options=None)`。
- **输入：** 正整数 contig 长度列表；可选 `genome_size` 以附加 NG50。
- **输出：** 计数、总长、最小/最大/平均长度、N50/L50、N90/L90 和排序后的长度。
- **示例代码：** `dnakit.assembly.calculate_stats([100, 60, 40], options={"genome_size": 250})`。
- **示例结果：** `{"result": {"contig_count": 3, "n50": 100, "ng50": 60}}`。

