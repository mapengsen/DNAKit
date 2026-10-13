# 原生 SAM 记录与比对统计

这些函数处理 SAM 文本或已经物化的 SAM-like 记录。它们不读取二进制 BAM/CRAM，不构建参考索引，也不替代 BWA、Bowtie2、minimap2 等外部比对器。比对记录统计方法可引用 SAMtools/BCFtools/HTSlib：<https://doi.org/10.1093/gigascience/giab008>。

## 1) `dnakit.alignment.parse_sam_records()`

- **作用：** 将 SAM 文本解析为可序列化的记录列表。
- **处理方法：** 读取表头和至少 11 列的 SAM 记录，解析 flag、坐标、CIGAR、序列、质量和标签。
- **API：** `dnakit.alignment.parse_sam_records(source)`。
- **输入：** SAM 文本或本地 SAM 文件路径。
- **输出：** 每条记录的字典列表；坐标转换为 0-based 起点。
- **示例代码：** `dnakit.alignment.parse_sam_records("r1\\t0\\tchr1\\t1\\t60\\t4M\\t*\\t0\\t0\\tATGC\\tIIII")`。
- **示例结果：** `[{"query_name": "r1", "reference": "chr1", "position": 0, "cigar": "4M"}]`。

## 2) `dnakit.alignment.calculate_coverage()`

- **作用：** 计算 SAM-like 记录的覆盖碱基数、总深度和平均深度。
- **处理方法：** 按 CIGAR 的 `M`、`=`、`X` 和参考消耗操作累计每个参考位置的深度。
- **API：** `dnakit.alignment.calculate_coverage(input)`。
- **输入：** `{"sam": "..."}` 或 `{"records": [...]}`。
- **输出：** `covered_bases`、`total_depth`、`mean_depth` 和分参考深度。
- **示例代码：** `dnakit.alignment.calculate_coverage({"sam": sam_text})`。
- **示例结果：** `{"covered_bases": 4, "total_depth": 4, "mean_depth": 1.0}`。

## 3) `dnakit.alignment.calculate_mapping_rate()`

- **作用：** 计算记录中已映射 reads 的比例。
- **处理方法：** 根据 SAM unmapped flag 或 `*` reference 统计。
- **API：** `dnakit.alignment.calculate_mapping_rate(input)`。
- **输入：** SAM 文本、路径或物化记录。
- **输出：** 0 到 1 之间的浮点比例。
- **示例代码：** `dnakit.alignment.calculate_mapping_rate({"records": records})`。
- **示例结果：** `0.75`。

## 4) `dnakit.alignment.generate_pileup()`

- **作用：** 生成按参考和位置排列的深度记录。
- **处理方法：** 复用本地 CIGAR 深度计算，不调用 samtools。
- **API：** `dnakit.alignment.generate_pileup(input)`。
- **输入：** SAM-like 记录。
- **输出：** `[{"reference": ..., "position": ..., "depth": ...}]`。
- **示例代码：** `dnakit.alignment.generate_pileup({"sam": sam_text})`。
- **示例结果：** `[{"reference": "chr1", "position": 0, "depth": 2}]`。

## 5) `dnakit.alignment.mark_duplicates()`

- **作用：** 按参考、位置、方向和序列标记重复记录。
- **处理方法：** 第一次出现保留，后续相同键的记录增加 `DUP=1` 标签。
- **API：** `dnakit.alignment.mark_duplicates(input)`。
- **输入：** SAM-like 记录。
- **输出：** 带标签的记录列表。
- **示例代码：** `dnakit.alignment.mark_duplicates({"records": records})`。
- **示例结果：** `[{"tags": {"DUP": "1"}}]`。

