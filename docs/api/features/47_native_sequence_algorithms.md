# 原生序列 k-mer、motif 与复杂度算法

这些函数只处理显式 DNA 字符串，不读取参考数据库，也不加载训练模型。方法来源包括 SeqKit2 的序列工具论文 <https://doi.org/10.1002/imt2.191>、MEME Suite 的 motif 分析论文 <https://doi.org/10.1093/nar/gkp335> 和 VSEARCH 的序列聚类论文 <https://doi.org/10.7717/peerj.2584>。

## 1) `dnakit.sequence.count_kmers()`

- **作用：** 统计序列中每种 k-mer 的出现次数。
- **处理方法：** 按重叠窗口扫描规范化 DNA，保留输入顺序无关的计数映射。
- **API：** `dnakit.sequence.count_kmers(sequence, k=5)`。
- **输入：** DNA 字符串和正整数 `k`。
- **输出：** k-mer 到计数的字典。
- **示例代码：** `dnakit.sequence.count_kmers("ATAT", k=2)`。
- **示例结果：** `{"AT": 2, "TA": 1}`。

## 2) `dnakit.sequence.count_canonical_kmers()`

- **作用：** 合并正向和反向互补 k-mer 的计数。
- **处理方法：** 对每个窗口取字典序较小的正向/反向互补表示。
- **API：** `dnakit.sequence.count_canonical_kmers(sequence, k=5)`。
- **输入：** DNA 字符串和正整数 `k`。
- **输出：** canonical k-mer 到计数的字典。
- **示例代码：** `dnakit.sequence.count_canonical_kmers("ATGC", k=2)`。
- **示例结果：** `{"AT": 1, "CA": 1, "GC": 1}`。

## 3) `dnakit.sequence.build_minimizer_index()`

- **作用：** 为序列构建有界 minimizer 索引。
- **处理方法：** 在 k-mer 窗口内选择字典序最小的 k-mer，并记录起点位置。
- **API：** `dnakit.sequence.build_minimizer_index(sequence, k=15, window=10)`。
- **输入：** DNA 字符串、k-mer 长度和 minimizer 窗口大小。
- **输出：** minimizer 到位置列表的字典。
- **示例代码：** `dnakit.sequence.build_minimizer_index("ACGTACGT", k=3, window=2)`。
- **示例结果：** `{"ACG": [0, 4]}`。

## 4) `dnakit.sequence.generate_kmer_table()`

- **作用：** 为多条序列生成带序列索引的 k-mer 表。
- **处理方法：** 逐条调用确定性的 k-mer 计数逻辑并展开为记录。
- **API：** `dnakit.sequence.generate_kmer_table(sequences, k=5)`。
- **输入：** DNA 字符串列表和 `k`。
- **输出：** 含 `sequence_index`、`kmer`、`count` 的记录列表。
- **示例代码：** `dnakit.sequence.generate_kmer_table(["ATAT"], k=2)`。
- **示例结果：** `[{"sequence_index": 0, "kmer": "AT", "count": 2}]`。

## 5) `dnakit.sequence.scan_pwm()`

- **作用：** 扫描序列并返回达到阈值的 PWM 命中。
- **处理方法：** 对每个窗口按 PWM 列分数求和，再与 `threshold` 比较。
- **API：** `dnakit.sequence.scan_pwm(sequence, pwm, threshold=0.0)`。
- **输入：** DNA 字符串、按位置排列的碱基分数映射和阈值。
- **输出：** 包含起点、窗口序列和分数的命中列表。
- **示例代码：** `dnakit.sequence.scan_pwm("ATGC", [{"A": 1}, {"T": 1}], threshold=1.5)`。
- **示例结果：** `[{"start": 0, "end": 2, "sequence": "AT", "score": 2.0}]`。

## 6) `dnakit.sequence.mask_low_complexity()`

- **作用：** 将低熵窗口替换为指定 mask 字符。
- **处理方法：** 计算滑动窗口 Shannon entropy，低于阈值的窗口替换为 `N`。
- **API：** `dnakit.sequence.mask_low_complexity(sequence, window=10, entropy_threshold=1.2)`。
- **输入：** DNA 字符串、窗口大小和熵阈值。
- **输出：** mask 后的 DNA 字符串。
- **示例代码：** `dnakit.sequence.mask_low_complexity("AAAAAAAA", window=4, entropy_threshold=1.2)`。
- **示例结果：** `"NNNNNNNN"`。

## 7) `dnakit.sequence.cluster_sequences_by_identity()`

- **作用：** 按确定性的序列身份阈值生成聚类。
- **处理方法：** 使用等长位置身份比例，把新序列归入第一个达到阈值的代表簇。
- **API：** `dnakit.sequence.cluster_sequences_by_identity(sequences, threshold=0.9)`。
- **输入：** DNA 字符串列表和 0 到 1 的身份阈值。
- **输出：** 每个簇的输入序列索引列表。
- **示例代码：** `dnakit.sequence.cluster_sequences_by_identity(["ACGT", "ACGA"], threshold=0.75)`。
- **示例结果：** `[[0, 1]]`。

