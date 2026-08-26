# API 简化审计

## 结论

普通用户只需记住以下规则：

```python
import dnakit
from dnakit.ops import delete

dna = dnakit.DNA("ACGT", id="seq-1")
many = dnakit.DNA(["ACGT", "TTAA"])
loaded = dnakit.read("input.fa", mode="dna")
report = dnakit.validate(dna)
edited = delete(dna, 1, 2)
```

单条、多条、下标、切片和常用编辑都保持 `DNA` 类型。内部精确类型和详细审计结果仍存在，但不要求普通用户先学习。

## 已统一的重复入口

| 场景 | 普通用户入口 | 保留的高级/兼容入口 | 处理结果 |
|---|---|---|---|
| 构造序列、记录、数据集 | `DNA(...)` | `DNASequence`、`DNARecord`、`DNASet`、`from_sequences()`、`from_records()` | 单条和多条始终返回 `DNA` |
| 从多条中选择 | `dna[index]`、`dna[slice]` | `dna.records`、`dna.record`、`dna.dataset` | 用户选择后仍得到 `DNA` |
| 文件读写 | `read(..., mode="dna"|"stream")`、`write(...)` | `read_one()`、`read_set()`、显式索引 API | 用参数选择物化对象或大文件流；`write()` 可直接消费流，IO-005 仅保留为需求追踪编号 |
| 合法性检查 | `validate(...)` | `validate_set()` | 根据 `DNA` 内记录数自动返回单条或集合报告 |
| 序列/记录编辑 | 无后缀的 `insert/delete/substitute/mask/trim` | 对应 `*_record()` | 输入 `DNA` 时同步注释并返回 `DNA` |
| 方向和环状操作 | `reverse_complement/rotate/canonical_origin` | 对应 `*_record()` | 输入 `DNA` 时同步注释并返回 `DNA` |
| 元数据操作 | 原有 `with_metadata/merge_metadata/filter_by_metadata/select_metadata` | 无需新名称 | 直接接受并返回 `DNA` |
| 常用分析 | 原有描述符、指纹、模式、相似度、比对、评价和可视化名称 | 内部序列/记录输入仍兼容 | 直接接受 `DNA`，不要求手动拆对象 |

## 没有合并的名称

以下名称虽然相关，但语义不同，因此不能只靠输入类型猜测：

- `concat()` 与 `concat_overlap()`：普通拼接与去重叠拼接不是同一算法。
- 精确去重、近似去重和基于模型的聚类：算法、阈值和资源边界不同。
- FASTA/FASTQ 记录 I/O 与 GFF3/BED/AGP document codec：数据结构和格式规则不同。
- `reverse()`、`complement()`、`reverse_complement()`：代表三个不同的方向变换。
- 各领域分析函数：结果含义不同，保留明确的科学名称，避免“万能 analyze()”隐藏算法。

## 兼容边界

本次没有删除或重命名旧公开 API。旧代码可继续使用明确分层的核心对象和 `*_record()` 审计接口；新文档只把统一入口放到主学习路径。
